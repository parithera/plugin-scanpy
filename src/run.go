package js

import (
	"log"
	"os"
	"os/exec"
	"path"
	"path/filepath"
	"time"

	codeclarity "github.com/CodeClarityCE/utility-types/codeclarity_db"
	exceptionManager "github.com/CodeClarityCE/utility-types/exceptions"
	"github.com/uptrace/bun"

	"github.com/parithera/plugin-scanpy/src/types"
	"github.com/parithera/plugin-scanpy/src/utils/output_generator"
)

// Start is a function that analyzes the source code directory and generates a software bill of materials (SBOM) output.
// It returns an sbomTypes.Output struct containing the analysis results.
func Start(sourceCodeDir string, codeclarityDB *bun.DB) types.Output {
	return ExecuteScript(sourceCodeDir)
}

func ExecuteScript(sourceCodeDir string) types.Output {
	start := time.Now()

	files, err := filepath.Glob(sourceCodeDir + "/*.fastq.gz")
	if err != nil {
		log.Fatal(err)
	}

	h5_files, err := filepath.Glob(sourceCodeDir + "/*.h5")
	if err != nil {
		log.Fatal(err)
	}

	if len(files) == 0 && len(h5_files) == 0 {
		return generate_output(start, "no fastq file", codeclarity.SUCCESS, []exceptionManager.Error{})
	}

	outputPath := path.Join(sourceCodeDir, "scanpy")
	os.MkdirAll(outputPath, os.ModePerm)

	starPath := path.Join(sourceCodeDir, "STAR", "outSolo.out", "Gene", "filtered")
	args := []string{"scripts/main.py", starPath}

	// If there is no fastq then it means there is an h5
	if len(files) == 0 {
		starPath = path.Join(sourceCodeDir, "data.h5")
		args = []string{"scripts/main_h5.py", starPath}
	}

	// Run Rscript in sourceCodeDir
	cmd := exec.Command("python3", args...)
	_, err = cmd.CombinedOutput()
	if err != nil {
		// panic(fmt.Sprintf("Failed to run Rscript: %s", err.Error()))
		codeclarity_error := exceptionManager.Error{
			Private: exceptionManager.ErrorContent{
				Description: err.Error(),
				Type:        exceptionManager.GENERIC_ERROR,
			},
			Public: exceptionManager.ErrorContent{
				Description: "The script failed to execute",
				Type:        exceptionManager.GENERIC_ERROR,
			},
		}
		return generate_output(start, nil, codeclarity.FAILURE, []exceptionManager.Error{codeclarity_error})
	}

	return generate_output(start, "done", codeclarity.SUCCESS, []exceptionManager.Error{})
}

func generate_output(start time.Time, data any, status codeclarity.AnalysisStatus, errors []exceptionManager.Error) types.Output {
	formattedStart, formattedEnd, delta := output_generator.GetAnalysisTiming(start)

	output := types.Output{
		Result: types.Result{
			Data: data,
		},
		AnalysisInfo: types.AnalysisInfo{
			Errors: errors,
			Time: types.Time{
				AnalysisStartTime: formattedStart,
				AnalysisEndTime:   formattedEnd,
				AnalysisDeltaTime: delta,
			},
			Status: status,
		},
	}
	return output
}
