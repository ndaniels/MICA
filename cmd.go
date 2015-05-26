package cablastp

import (
	"bytes"
	"fmt"
	"os/exec"
	"strings"
)

// Exec runs a command created with 'Command' in the os/exec package, and
// converts anything reported to stderr to a Go error value.
//
// Note that if the command returns successfully, the error is guaranteed to
// be nil.
func Exec(cmd *exec.Cmd) error {
	var stderr bytes.Buffer

	cmd.Stderr = &stderr
	fullCmd := strings.Join(cmd.Args, " ")

	Vprintf("%s\n", fullCmd)
	if err := cmd.Run(); err != nil {
		if stderr.Len() > 0 {
			return fmt.Errorf(
				"Error running '%s': '%s'. \n\nstderr:\n%s",
				fullCmd, err, stderr.String())
		}
		return fmt.Errorf("Error running '%s': '%s'.", fullCmd, err)
	}
	return nil
}
