package scanpy

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

// Start initiates the scanpy execution process.
// It takes the source code directory and a database connection as input.
// It returns the output of the scanpy execution.
func Start(sourceCodeDir string, codeclarityDB *bun.DB) types.Output {
	return ExecuteScript(sourceCodeDir)
}

func ExecuteScript(sourceCodeDir string) types.Output {
	startTime := time.Now()

	// Find fastq.gz files in the source directory.
	files, err := filepath.Glob(sourceCodeDir + "/*.fastq.gz")
	if err != nil {
		log.Fatal(err)
	}

	// Find h5 files in the source directory.
	h5_files, err := filepath.Glob(sourceCodeDir + "/*.h5")
	if err != nil {
		log.Fatal(err)
	}

	// If no fastq or h5 files are found, return an output indicating no files were found.
	if len(files) == 0 && len(h5_files) == 0 {
		return generate_output(startTime, "no fastq file", codeclarity.SUCCESS, []exceptionManager.Error{})
	}

	// Create the output directory.
	outputPath := path.Join(sourceCodeDir, "scanpy")
	os.MkdirAll(outputPath, os.ModePerm)

	// Define the default path for STAR output.
	starPath := path.Join(sourceCodeDir, "STAR", "outSolo.out", "Gene", "filtered")
	// Define the arguments for the script.
	args := []string{"scripts/main.py", starPath}

	// If no fastq files are found, assume h5 files are present and adjust the script and path accordingly.
	if len(files) == 0 {
		starPath = path.Join(sourceCodeDir, "data.h5")
		args = []string{"scripts/main_h5.py", starPath}
	}

	// Execute the python script.
	cmd := exec.Command("python3", args...)
	_, err = cmd.CombinedOutput()
	if err != nil {
		// Handle script execution errors.
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
		return generate_output(startTime, nil, codeclarity.FAILURE, []exceptionManager.Error{codeclarity_error})
	}

	return generate_output(startTime, "done", codeclarity.SUCCESS, []exceptionManager.Error{})
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
