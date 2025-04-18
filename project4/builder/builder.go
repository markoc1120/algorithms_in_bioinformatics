package builder

import (
	"fmt"
	"math"
	"strconv"
	"strings"

	"github.com/markoc1120/nj/utils"
)

type Node struct {
	ID       string
	Children map[*Node]float64
}

type Relation struct {
	Parent   string
	Distance float64
}

func (n *Node) AddChild(child *Node, distance float64) {
	n.Children[child] = distance
}

func NewNode(id string) *Node {
	return &Node{
		ID:       id,
		Children: make(map[*Node]float64),
	}
}

func sum(l []float64) float64 {
	sum := 0.0
	for _, el := range l {
		sum += el
	}
	return sum
}

func BuildNodes(root *Node, distanceMap map[string]Relation) *Node {
	nodeMap := make(map[string]*Node)
	nodeMap[root.ID] = root

	queue := []*Node{root}
	for len(queue) > 0 {
		var nodeList []*Node

		for _, node := range queue {
			for childID, rel := range distanceMap {
				if rel.Parent == node.ID {
					childNode := NewNode(childID)
					node.AddChild(childNode, rel.Distance)
					nodeMap[childID] = childNode
					nodeList = append(nodeList, childNode)
					delete(distanceMap, childID)
				}
			}
		}

		queue = nodeList
	}
	return root
}

func GetClosesstPair(
	n, activeCount int, active []bool, rowSums []float64, matrix [][]float64,
) (minI, minJ int) {
	minI, minJ, minQ := 0, 0, math.Inf(1)

	for i := range n {
		if !active[i] {
			continue
		}

		for j := i + 1; j < n; j++ {
			if !active[j] {
				continue
			}

			currQ := float64(activeCount-2)*matrix[i][j] - rowSums[i] - rowSums[j]

			if currQ < minQ {
				minQ, minI, minJ = currQ, i, j
			}
		}
	}
	return minI, minJ
}

func CalculateClosestPairEdges(
	minI, minJ, activeCount int, rowSums []float64, matrix [][]float64,
) (edgeI, edgeJ float64) {
	dIJ := matrix[minI][minJ]
	rI := rowSums[minI] / float64(activeCount-2)
	rJ := rowSums[minJ] / float64(activeCount-2)
	edgeI = 0.5 * (dIJ + rI - rJ)
	edgeJ = dIJ - edgeI
	return edgeI, edgeJ
}

func UpdateDistances(minI, minJ, n int, active []bool, rowSums []float64, matrix [][]float64) {
	rowSums[minI] = 0.0
	for k := range n {
		if !active[k] || k == minI || k == minJ {
			continue
		}

		newDist := 0.5 * (matrix[minI][k] + matrix[minJ][k] - matrix[minI][minJ])
		matrix[minI][k] = newDist
		matrix[k][minI] = newDist

		rowSums[minI] += newDist
		rowSums[k] = rowSums[k] + newDist - matrix[k][minI] - matrix[k][minJ]
	}
	active[minJ] = false
}

// TODO: test
func Termination(
	n int, active []bool, nodeIDs []string, matrix [][]float64, distanceMap map[string]Relation,
) *Node {
	var finalNodes []int
	for i := range n {
		if active[i] {
			finalNodes = append(finalNodes, i)
		}
	}
	a, b, c := finalNodes[0], finalNodes[1], finalNodes[2]
	dAB := matrix[a][b]
	dAC := matrix[a][c]
	dBC := matrix[b][c]

	edgeA := 0.5 * (dAB + dAC - dBC)
	edgeB := 0.5 * (dAB + dBC - dAC)
	edgeC := 0.5 * (dAC + dBC - dAB)
	rootID := nodeIDs[a] + nodeIDs[b] + nodeIDs[c]
	root := NewNode(rootID)

	distanceMap[nodeIDs[a]] = Relation{Parent: rootID, Distance: edgeA}
	distanceMap[nodeIDs[b]] = Relation{Parent: rootID, Distance: edgeB}
	distanceMap[nodeIDs[c]] = Relation{Parent: rootID, Distance: edgeC}
	return root
}

func NeiSaitou(phylip *utils.Phylip) *Node {
	var root *Node
	distanceMap := map[string]Relation{}
	n := len(phylip.IDs)

	// Base case
	if n <= 2 {
		switch {
		case n == 2:
			nodeI, nodeJ := phylip.IDs[0], phylip.IDs[1]
			distanceMap[nodeI] = Relation{Parent: nodeJ, Distance: phylip.Matrix[0][1]}
			root = NewNode(nodeJ)
		case n < 2:
			root = NewNode(phylip.IDs[0])
		}
		return BuildNodes(root, distanceMap)
	}

	matrix := phylip.Matrix
	active := make([]bool, n)
	for i := range active {
		active[i] = true
	}

	rowSums := make([]float64, n)
	for i := range rowSums {
		rowSums[i] = sum(matrix[i])
	}

	nodeIDs := make([]string, n)
	copy(nodeIDs, phylip.IDs)

	activeCount := n
	for activeCount > 3 {
		minI, minJ := GetClosesstPair(n, activeCount, active, rowSums, matrix)
		edgeI, edgeJ := CalculateClosestPairEdges(minI, minJ, activeCount, rowSums, matrix)

		// update distanceMap with the pairs
		node1, node2 := nodeIDs[minI], nodeIDs[minJ]
		newNodeID := node1 + node2
		distanceMap[node1] = Relation{Parent: newNodeID, Distance: edgeI}
		distanceMap[node2] = Relation{Parent: newNodeID, Distance: edgeJ}
		nodeIDs[minI] = newNodeID

		UpdateDistances(minI, minJ, n, active, rowSums, matrix)
		activeCount--
	}

	root = Termination(n, active, nodeIDs, matrix, distanceMap)
	return BuildNodes(root, distanceMap)
}

// TODO: test
func WriteNewickFormat(ids []string, root *Node) string {
	var postorderTraversal func(*Node) string
	postorderTraversal = func(node *Node) string {
		if len(node.Children) == 0 {
			// Check if this is a leaf node (original taxa)
			id, err := strconv.Atoi(node.ID)
			if err == nil && id >= 1 && id <= len(ids) {
				return ids[id-1]
			}
			return node.ID
		}

		// Process children
		var parts []string
		for child, distance := range node.Children {
			childStr := postorderTraversal(child)
			parts = append(parts, fmt.Sprintf("%s:%g", childStr, distance))
		}

		return "(" + strings.Join(parts, ",") + ")"
	}

	// Handle the case when root has 3 children
	if len(root.Children) == 3 {
		var children []*Node
		var distances []float64

		for child, distance := range root.Children {
			children = append(children, child)
			distances = append(distances, distance)
		}

		// Using all three children directly without any special formatting
		return "(" +
			postorderTraversal(children[0]) + ":" + fmt.Sprintf("%g", distances[0]) + "," +
			postorderTraversal(children[1]) + ":" + fmt.Sprintf("%g", distances[1]) + "," +
			postorderTraversal(children[2]) + ":" + fmt.Sprintf("%g", distances[2]) +
			");"
	}

	// Handle normal case
	return postorderTraversal(root) + ";"
}
