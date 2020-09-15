function [A,varargout] = WienerFilter_definitive(ToPredict,Predictors,varargin) 
% This Function take as input whatever predictors and data to be predicted
% and train a Wiener Filter with them
% Inputs are:
% - ToPredict = NxM matrix rows are observations, columns different element
% (can be cell array in case we have several good / bad junks interval)
% - ToBePredicted = NxM matrix rows are obs, columns features
% (can be cell array in case we have several good / bad junks interval)
% - SampleRate = scalar with sampling rate of the variables
% - TimePerPrediction = time interval between every prediction in mseconds

% - FeatureLength = amount in sec of previous time to be used for
%                   prediction
% - cross_val = values from 0 to 1, 0 no crossval otherwise the it's a
% percentage of data used for testing set
% - FeatureSelection = scalar, if 0 use all the features, whereas use a
%                      selection
% Outputs are:  A,Test_perf,MSE,pol
% - The matrix for the linear correlation Y = A * X
% - The performance on the 10% of the data left out as testing set:
% R values
% MSE values 
% in case of non crossvalidated results, the algorithm will be trained on
% 90% of the data and tested on the other dataset
%
% Posiibility to run cross-validation for the results
% Written by Nicolo' Macellari, Ecole Polytechnique de Lausanne
% nicolo.macellari@epfl.ch
% July 2020

% assignement of the variables.
if nargin < 3
    disp('missing sample Frequency, should i guess it?')
    return
else
    cross_val = 0;
    FeatureSelection = 0;
    FeatureLength = 0.5; 
    new_samp = 50;
    switch nargin
        case 3
            SampleRate = varargin{1};          
            disp('newSampFq 50, feature length 0.5')
        case 4
            SampleRate = varargin{1};
            new_samp = varargin{2};            
            disp('feature length 0.5')
        case 5
            SampleRate = varargin{1};
            new_samp = varargin{2};
            FeatureLength = varargin{3};
        case 6
            SampleRate = varargin{1};
            new_samp = varargin{2};
            FeatureLength = varargin{3};
            cross_val = varargin{4};
        case 7
            SampleRate = varargin{1};
            new_samp = varargin{2};
            FeatureLength = varargin{3};
            cross_val = varargin{4};
    end
end

resemp = SampleRate / new_samp;
f_plot = 1;


%% Check for junks
if iscell(ToPredict)
    
    X = [];
    Y = [];
    
    % for loop for each trunk
    for ii = 1:size(ToPredict,2)
        if isempty(ToPredict{1,ii})
            continue
        end
        Bin_EMG = [];
        Bin_spike = [];
        X_cut = [];
        Y_cut = [];
        
        % checks length
        Data_toPredict = ToPredict{1,ii}(1:length(Predictors{1,ii}),:);
        Data_Predictors = Predictors{1,ii};
          
        % resampling to right samp fq
        x = 0 : 1 : length(Data_toPredict)-1;
        xq = 0 : resemp : length(Data_toPredict)-1;
        Bin_EMG = interp1(x,Data_toPredict,xq);
        Bin_spike = interp1(x,Data_Predictors,xq);
        
        % Preparing X and Y for training
        SamplePredictionRate = FeatureLength  *  new_samp;
        
        j = SamplePredictionRate;
        a = 1;
        while j < size(Bin_spike,1)
            x = Bin_spike(a:j,:)';
            X_cut(a,:) = x(:)';
            Y_cut(a,:) = Bin_EMG(j+1,:);
            j = j+1;
            a = a+1;
        end
        
         if ~isempty (X_cut)
%             [XX,XmatMean,XmatStd] = sphereData(X_cut');
%             YY = sphereData(Y_cut');
%             X_cut = XX{1,1}';
%             Y_cut = YY{1,1}';
%             
             X = [X; X_cut];
             Y = [Y;Y_cut];
         end
        
    end
    
else
    
    
    % Feature selection
    if FeatureSelection ~= 0
        Gr = zeros(size(ToPredict,2),size(ToPredict,1));
        for i = 1: size(ToPredict,2)
            Gr(i,:) = (ToPredict(:,i)/max(ToPredict(:,i)) > 0.2)';
        end
        
        for i = 1: size(ToPredict,2)
            inmodel(i,:) = rankfeatures(Predictors',Gr(i,:));
        end
    else
        index = 1:size(Predictors,2);
        inmodel = repmat(index,size(ToPredict,2),1);
    end
    
    
    %% resampling to right samp fq
    
    x = 0 : 1 : length(ToPredict)-1;
    xq = 0 : resemp : length(ToPredict)-1;
    Bin_EMG = interp1(x,ToPredict,xq);
    Bin_spike = interp1(x,Predictors,xq);
    
    
    %% Preparing X and Y for training
    SamplePredictionRate = FeatureLength  *  new_samp;
    
    j = SamplePredictionRate;
    a = 1;
    while j < size(Bin_spike,1)
        x = Bin_spike(a:j,:)';
        X(a,:) = x(:)';
        Y(a,:) = Bin_EMG(j+1,:);
        j = j+1;
        a = a+1;
    end
    
    
end

% for sess = length(Predictors)
%
% end

[XX,XmatMean,XmatStd] = sphereData(X);
YY = sphereData(Y);
X = XX{1,1};
Y = YY{1,1};

finalR2 = [];
% cross-validation of the result to have a more consistent evaluation
if cross_val ~=0
    ridge_regression = [0 power(10,linspace(-1,7,20))];

    for reg = 1:length(ridge_regression)
        lambda = ridge_regression(reg); % values for L2 regularization
        Cross_val_wiener_innerFunction;
        finalR2(reg,:) = Test_perf;
    end

    [x,idx] = max(finalR2(:,5));
    lambda = ridge_regression(idx);
    Full_Data_set;
    
    varargout{1} = finalR2(idx,:);
    varargout{2} = MSE;
    varargout{3} = pol;
    
    
    return
end


% no cross_val with plots
No_cross_val_training


varargout{1} = Test_perf;
varargout{2} = MSE;
varargout{3} = pol;
%varargout{4} = mu;
varargout{5} = VAF;

end

