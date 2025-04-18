package builder_test

import (
	"math"
	"testing"

	"github.com/markoc1120/nj/builder"
	"github.com/markoc1120/nj/utils"
)

type testCase struct {
	n, activeCount, expectedMinI, expectedMinJ int
	expectedEdgeI, expectedEdgeJ               float64
	active, expectedActive                     []bool
	rowSums, expectedRowSums                   []float64
	phylip                                     utils.Phylip
	expectedMatrix                             [][]float64
}

func TestBuilder(t *testing.T) {
	testCase := testCase{
		n:               3,
		activeCount:     3,
		expectedMinI:    0,
		expectedMinJ:    1,
		expectedEdgeI:   0.0,
		expectedEdgeJ:   0.11,
		active:          []bool{true, true, true},
		expectedActive:  []bool{true, false, true},
		rowSums:         []float64{0.33, 0.44, 0.55},
		expectedRowSums: []float64{0.22, 0.44, 0.22},
		phylip: utils.Phylip{
			IDs: []string{"A", "B", "C"},
			Matrix: [][]float64{
				{0.00, 0.11, 0.22},
				{0.11, 0.00, 0.33},
				{0.22, 0.33, 0.00},
			},
		},
		expectedMatrix: [][]float64{
			{0, 0.11, 0.22},
			{0.11, 0, 0.33},
			{0.22, 0.33, 0},
		},
	}

	t.Run("test NeiSaitou algorithm", func(t *testing.T) {
		got := builder.NeiSaitou(&testCase.phylip)
		if got.ID != "AB" && got.ID != "ABC" && got.ID != "CA" {
			t.Errorf("Expected root ID to be an internal node, got %s", got.ID)
		}

		nodeCount := countNodes(got)
		expectedNodeCount := 4
		if nodeCount != expectedNodeCount {
			t.Errorf("Expected %d nodes in tree, got %d", expectedNodeCount, nodeCount)
		}

		foundTaxa := findAllLeafIDs(got)
		expectedTaxa := map[string]bool{"A": true, "B": true, "C": true}
		for id := range expectedTaxa {
			if !foundTaxa[id] {
				t.Errorf("Expected to find taxon %s in tree, but it's missing", id)
			}
		}
	})

	t.Run("test GetClosesstPair function", func(t *testing.T) {
		minI, minJ := builder.GetClosesstPair(
			testCase.n,
			testCase.activeCount,
			testCase.active,
			testCase.rowSums,
			testCase.phylip.Matrix,
		)
		if testCase.expectedMinI != minI || testCase.expectedMinJ != minJ {
			t.Errorf(
				"got minI: %d minJ: %d, want minI: %d minJ: %d",
				minI,
				minJ,
				testCase.expectedMinI,
				testCase.expectedMinJ,
			)
		}
	})

	t.Run("test CalculateClosestPairEdges function", func(t *testing.T) {
		edgeI, edgeJ := builder.CalculateClosestPairEdges(
			testCase.expectedMinI,
			testCase.expectedMinJ,
			testCase.activeCount,
			testCase.rowSums,
			testCase.phylip.Matrix,
		)
		if !roughlyEqualFloat64(edgeI, testCase.expectedEdgeI) ||
			!roughlyEqualFloat64(edgeJ, testCase.expectedEdgeJ) {
			t.Errorf(
				"got edgeI: %f edgeJ: %f, want edgeI: %f edgeJ: %f ",
				edgeI,
				edgeJ,
				testCase.expectedEdgeI,
				testCase.expectedEdgeJ,
			)
		}
	})

	t.Run("test UpdateDistances function", func(t *testing.T) {
		builder.UpdateDistances(
			testCase.expectedMinI,
			testCase.expectedMinJ,
			testCase.n,
			testCase.active,
			testCase.rowSums,
			testCase.phylip.Matrix,
		)

		for i := range testCase.n {
			if !roughlyEqualFloat64(testCase.rowSums[i], testCase.expectedRowSums[i]) ||
				testCase.active[i] != testCase.expectedActive[i] {
				t.Errorf(
					"got rowSums: %+v active: %+v, want rowSums: %+v active: %+v",
					testCase.rowSums,
					testCase.active,
					testCase.expectedRowSums,
					testCase.expectedActive,
				)
			}
			for j := range testCase.n {
				if !roughlyEqualFloat64(
					testCase.phylip.Matrix[i][j],
					testCase.expectedMatrix[i][j],
				) {
					t.Errorf(
						"got matrix: %+v , want matrix: %+v",
						testCase.phylip.Matrix,
						testCase.expectedMatrix,
					)
				}
			}
		}
	})
}

func roughlyEqualFloat64(a, b float64) bool {
	const equalityThreshold = 1e-7
	return math.Abs(a-b) < equalityThreshold
}

func countNodes(node *builder.Node) int {
	if node == nil {
		return 0
	}

	count := 1
	for child := range node.Children {
		count += countNodes(child)
	}
	return count
}

func findAllLeafIDs(node *builder.Node) map[string]bool {
	result := make(map[string]bool)
	var traverse func(*builder.Node)

	traverse = func(n *builder.Node) {
		if n == nil {
			return
		}

		if len(n.Children) == 0 {
			result[n.ID] = true
		}

		for child := range n.Children {
			traverse(child)
		}
	}

	traverse(node)
	return result
}
