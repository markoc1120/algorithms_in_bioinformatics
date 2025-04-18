package utils_test

import (
	"reflect"
	"testing"
	"testing/fstest"

	"github.com/markoc1120/nj/utils"
)

const testPhyContent = `3
A 0.00 0.11 0.22
B 0.11 0.00 0.33
C 0.22 0.33 0.00`

func TestPhylipParse(t *testing.T) {
	t.Run("successful parse", func(t *testing.T) {
		fileSystem := fstest.MapFS{
			"test.phy": {Data: []byte(testPhyContent)},
		}

		got, _ := utils.ParsePhylip(fileSystem)
		expected := utils.Phylip{
			IDs: []string{"A", "B", "C"},
			Matrix: [][]float64{
				{0.00, 0.11, 0.22},
				{0.11, 0.00, 0.33},
				{0.22, 0.33, 0.00},
			},
			Filename: "test.phy",
		}

		if !reflect.DeepEqual(got[0], expected) {
			t.Errorf("got %+v, but expected %+v", got[0], expected)
		}
	})
}
