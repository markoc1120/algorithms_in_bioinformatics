package main

import (
	"encoding/csv"
	"fmt"
	"log"
	"os"
	"os/exec"
	"strings"
	"time"

	"github.com/markoc1120/nj/builder"
	"github.com/markoc1120/nj/utils"
)

type result struct {
	name     string
	duration string
	method   string
}

func main() {
	phylips, err := utils.ParsePhylip(os.DirFS("data"))
	if err != nil {
		fmt.Printf("Error parsing PHYLIP files in 'data' folder: %v\n", err)
		os.Exit(1)
	}

	quickTreePath := "./quicktree"
	resultPath := "result"
	timeLogs := [][]string{}
	resultChannel := make(chan result)

	if _, err := os.Stat(resultPath); os.IsNotExist(err) {
		err := os.Mkdir(resultPath, 0777)
		if err != nil {
			fmt.Printf("Couldn't create 'result' folder %v\n", err)
			os.Exit(1)
		}
	}
	for _, phylip := range phylips {
		go func(p utils.Phylip) {
			fmt.Printf("Running NeiSaitou: %s\n", p.Filename)

			start := time.Now()
			root := builder.NeiSaitou(&p)
			newick := builder.WriteNewickFormat(p.IDs, root)
			elapsed := time.Since(start)

			resultChannel <- result{p.Filename, fmt.Sprintf("%.2f", elapsed.Seconds()), "own_implementation"}

			err = os.WriteFile(
				fmt.Sprintf("result/ns_%s.txt", strings.TrimSuffix(p.Filename, ".phy")), []byte(newick), 0644,
			)
			if err != nil {
				fmt.Printf("Error writing Newick file: %v\n", err)
				os.Exit(1)
			}
			fmt.Printf("Finished NeiSaitou: %s\n", p.Filename)
		}(phylip)

		go func(p utils.Phylip) {
			fmt.Printf("Running QuickTree: %s\n", p.Filename)

			start := time.Now()
			cmd := exec.Command(quickTreePath, "-in", "m", fmt.Sprintf("data/%s", p.Filename))
			err := cmd.Run()
			elapsed := time.Since(start)

			resultChannel <- result{p.Filename, fmt.Sprintf("%.2f", elapsed.Seconds()), "QuickTree"}

			if err != nil {
				fmt.Printf("Error running QuickTree for %s: %v\n", p.Filename, err)
				os.Exit(1)
			}
			fmt.Printf("Finished Quicktree: %s\n", p.Filename)
		}(phylip)
	}

	for range len(phylips) * 2 {
		r := <-resultChannel
		timeLogs = append(timeLogs, []string{r.name, r.duration, r.method})
	}

	f, err := os.Create("result/times.csv")
	if err != nil {
		log.Fatalln("error creating file:", err)
	}

	w := csv.NewWriter(f)
	w.WriteAll(timeLogs)

	if err := w.Error(); err != nil {
		log.Fatalln("error writing csv:", err)
	}
}
