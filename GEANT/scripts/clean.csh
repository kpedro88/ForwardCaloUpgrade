#!/bin/tcsh

mv *.condor logs/
mv *.log logs/
mv *.std* logs/
mv G4_*.jdl logs/

foreach f ( wlyso pbwo4 hcal whcal )
foreach s ( in sh )
foreach i (`ls ${f}*.${s} | grep -v "temp" | grep -v "test"` )
mv $i runs/
end
end
end
