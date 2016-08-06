function [out_struct] = elmMOD(regr_in,regr_out,NumberofHiddenNeurons,ActivationFunction)


    %%%%    Authors:    MR QIN-YU ZHU AND DR GUANG-BIN HUANG
    %%%%    NANYANG TECHNOLOGICAL UNIVERSITY, SINGAPORE
    %%%%    EMAIL:      EGBHUANG@NTU.EDU.SG; GBHUANG@IEEE.ORG
    %%%%    WEBSITE:    http://www.ntu.edu.sg/eee/icis/cv/egbhuang.htm
    %%%%    DATE:       APRIL 2004


    
% The original file was modified for the integration with the
% Timbre Space Analyzer and Mapper (TSAM).
% 
% The TSAM is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% The TSAM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Less General Public License for more details.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with TSAM.  If not, see <http://www.gnu.org/licenses/>.
%
% The TSAM can be obtained at http://stefanofasciani.com/tsam.html
% TSAM Copyright (C) 2016 Stefano Fasciani, University of Wollongong
% Inquiries: stefanofasciani@stefanofasciani.com

    
%%%%%%%%%%% Load training dataset
T=regr_out;
P=regr_in;
clear train_data;                                   %   Release raw training data array

NumberofTrainingData=size(P,2);
NumberofInputNeurons=size(P,1);                                              %   end if of Elm_Type

%%%%%%%%%%% Random generate input weights InputWeight (w_i) and biases BiasofHiddenNeurons (b_i) of hidden neurons
InputWeight=rand(NumberofHiddenNeurons,NumberofInputNeurons)*2-1;
BiasofHiddenNeurons=rand(NumberofHiddenNeurons,1);
tempH=InputWeight*P;
clear P;                                            %   Release input of training data 
ind=ones(1,NumberofTrainingData);
BiasMatrix=BiasofHiddenNeurons(:,ind);              %   Extend the bias matrix BiasofHiddenNeurons to match the demention of H
tempH=tempH+BiasMatrix;

%%%%%%%%%%% Calculate hidden neuron output matrix H
switch lower(ActivationFunction)
    case {'sig','sigmoid'}
        %%%%%%%% Sigmoid 
        H = 1 ./ (1 + exp(-tempH));
    case {'sin','sine'}
        %%%%%%%% Sine
        H = sin(tempH);    
    case {'hardlim'}
        %%%%%%%% Hard Limit
        %H_test = hardlim(tempH);
        H = double(tempH >= 0);
    case {'tribas'}
        %%%%%%%% Triangular basis function
        %H = tribas(tempH);
        H = max(0,1-abs(tempH));
    case {'radbas'}
        %%%%%%%% Radial basis function
        H = exp(-(tempH.*tempH));
        %%%%%%%% More activation functions can be added here                
end
clear tempH;                                        %   Release the temparary array for calculation of hidden neuron output matrix H

%%%%%%%%%%% Calculate output weights OutputWeight (beta_i)
OutputWeight=pinv(H') * T';                        % implementation without regularization factor //refer to 2006 Neurocomputing paper

%%%%%%%%%%% Calculate the training accuracy
Y=(H' * OutputWeight)';                             %   Y: the actual output of the training data

TrainingAccuracy=sqrt(mse(T - Y));               %   Calculate training accuracy (RMSE) for regression case

clear H;

out_struct.accuracy=TrainingAccuracy;
out_struct.outweight=OutputWeight;
out_struct.inweight=InputWeight;
out_struct.bias=BiasofHiddenNeurons;
out_struct.actfunct=ActivationFunction;


