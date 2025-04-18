package utils

import (
	"bufio"
	"io"
	"io/fs"
	"strconv"
	"strings"
)

type Phylip struct {
	IDs      []string
	Matrix   [][]float64
	Filename string
}

func ParsePhylip(fileSystem fs.FS) ([]Phylip, error) {
	dir, err := fs.ReadDir(fileSystem, ".")
	if err != nil {
		return nil, err
	}

	var phylips []Phylip

	for _, file := range dir {
		phylip, err := getPhylip(fileSystem, file.Name())
		if err != nil {
			return nil, err
		}
		phylips = append(phylips, phylip)
	}
	return phylips, nil
}

func getPhylip(fileSystem fs.FS, fileName string) (Phylip, error) {
	phylipFile, err := fileSystem.Open(fileName)
	if err != nil {
		return Phylip{}, err
	}
	defer phylipFile.Close()
	return readPhylip(phylipFile, fileName)
}

func readPhylip(phylipFile io.Reader, fileName string) (Phylip, error) {
	scanner := bufio.NewScanner(phylipFile)

	scanner.Scan()
	count, err := strconv.Atoi(scanner.Text())
	if err != nil {
		return Phylip{}, err
	}

	ids := []string{}
	matrix := [][]float64{}

	readLine := func() error {
		scanner.Scan()
		line := scanner.Text()
		fields := strings.Fields(line)
		ids = append(ids, fields[0])

		distances := []float64{}
		for _, val := range fields[1:] {
			distance, err := strconv.ParseFloat(val, 64)
			if err != nil {
				return err
			}
			distances = append(distances, distance)
		}
		matrix = append(matrix, distances)
		return nil
	}
	for range count {
		err := readLine()
		if err != nil {
			return Phylip{}, err
		}
	}

	return Phylip{
		IDs:      ids,
		Matrix:   matrix,
		Filename: fileName,
	}, nil
}
