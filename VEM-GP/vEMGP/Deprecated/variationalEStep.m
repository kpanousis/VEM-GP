function  [seq , varBound]=variationalEStep(seq,params,InferenceMethod)
%VARIATIONALESTEP File implementing the E Step
%Konstantinos Panagiotis Panousis
%Gatsby Computational Neuroscience Unit
%University College London
%8 June 2015
%Used to call the Inference Method of PopSpikeDyn package in order to
%perform inference
%Inputs: seq: Containing fields
%              seq.y: NxT for a spike raster from N neurons
%              seq.T: the 2nd dimension of seq.Y
%        params:  the parameters of the GLM
%        InferenceMethod: the inference method

 hypers=params.model.hypers;
 params.model.Q0=exp(hypers(1));
 params.model.B=1;
 params.model.A=exp((1/samplingRate)*hypers(4));
 
[seq, varBound]=InferenceMethod(params,seq);
end

