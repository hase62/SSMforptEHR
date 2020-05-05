# SSMforptEHR

## There exist two programs for the prediction of blood test values using state space models with Gaussian and skew-t distributed observation noises. 

### SSM with Gaussian noise requires the following arguments. 

<main: mainVARSSMEx.java>

arg1: Root directory

arg2: Observation data

arg3: Given regularoty structure among hidden variables

arg4: Id

arg5: Setting file

arg6: Output directory

arg7: External input

arg8: External input at t=0

arg9: Given regulatory structure from external input

## For example, you can run

java -jar SSM.jar testdata_1 PG.Obs.txt NO 0 varssm_test.set Output PG.Drug.txt Z0.txt NO

### SSM with skew-t distributed noise requires the following arguments. 

<main: mainstFilter.java>

arg1: Root directory

arg2: Observation data

arg3: Given regularoty structure among hidden variables

arg4: Id

arg5: Setting file

arg6: Output directory

arg7: External input

arg8: External input at t=0

arg9: Given regulatory structure from external input

arg10: Parameter nu for Gamma distribution

## For example, you can indicate

java -jar stSSM.jar testdata_1 PG.Obs.txt NO 0 varssm_test.set Output PG.Drug.txt Z0.txt NO 5

java -jar stSSM.jar testdata_2 Linear/obsData_norm_sk_pred_8nodes_Noise0.1-0.3ts1.0ID0.txt NO 0 varssm_linear.set obsData_norm_sk_pred_8nodes_Noise0.1-0.3ts1.0ID0 NO NO NO 5
