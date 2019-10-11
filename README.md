# SSMforptEHR

This program requires the following arguments. 

java -jar VARSSM.jar <root directory> <observation data> <given regularoty structure among hidden variables> <id> <setting file> <output directory> <external input> <external input at t=0> <given regulatory structure from external input>

For example, you can run it as ...

testdata PG.Obs.txt NO 0 varssm_test.set Output PG.Drug.txt Z0.txt NO 
