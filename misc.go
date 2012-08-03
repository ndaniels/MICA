package cablastp

import (
	"flag"
	"fmt"
)

func PrintFlagDefaults() {
	flag.VisitAll(func(fg *flag.Flag) {
		fmt.Printf("--%s=\"%s\"\n\t%s\n", fg.Name, fg.DefValue, fg.Usage)
	})
}
