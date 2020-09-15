% scripts for the cross-validation of the Wiener

length_segments = length(X(:,1)) * cross_val; % 10 percent
indeces = ceil(1:length_segments:length(X(:,1)));

% removing zero channels for matrix inversion and manipulaiton
num_zeros = sum(X,1);
X(:,num_zeros == 0) = [];

for j = 1:length(indeces)-1
    clear X_training Y_training X_test  A pol test training
    
    % finding indices
    test = indeces(j) : indeces(j+1);
    training = 1:max(indeces);
    training(test) = [];
    
    % division in training set and test set
    X_training = X(training,:);
    Y_training = Y(training,:);
    X_test = X(test,:);
    Y_test{j} = Y(test,:);
    
    % Matrix calculation
    X_training = [ones(size(X_training,1),1) X_training];
    X_test = [ones(size(X_test,1),1) X_test];
    Y_tr = Y_training;
    
    A = inv(transpose(X_training)*X_training + lambda *eye(size(X_training,2)))...
        *transpose(X_training)*Y_training;
    
    % Polynomial fitting (only one iteration)
    Y_rs = X_training * A;
    for ii = 1:size(Y_tr,2)
        [pol(ii,:),s] = polyfit(Y_rs(:,ii), Y_tr(:,ii),2);
    end
    
    Y_res = X_test * A;
    Y_fin{j} = [];
    
    % Polynomial application to the output of Wiener filter
    for i=1:size(Y_res,2)
        Y_fin{j}(:,i) = polyval(pol(i,:),Y_res(:,i));
    end
    
end

% concateneting segments
YfinPredicted = [];
YfinReal = [];
for j = 1:length(Y_fin)
    YfinReal = [YfinReal; Y_test{j}];
    YfinPredicted = [YfinPredicted; Y_fin{j}];
end

for k = 1:size(YfinPredicted,2)
        R = (corrcoef(YfinReal(:,k),YfinPredicted(:,k)));
        %Test_perf(k) = R(1,2);
        Test_perf(k) = 1 - var(YfinReal(:,k) - YfinPredicted(:,k))/var(YfinReal(:,k));
        MSE(k) = sum((YfinReal(:,k) - YfinPredicted(:,k)).^2) / size(YfinReal,1);
end

% clearing for full dataset training
clear pol Y_fin Y_res A mupol


