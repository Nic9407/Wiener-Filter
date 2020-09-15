% standard training no cross val, last 10% data taken as test set
% removing zero channels for matrix inversion and manipulaiton
num_zeros = sum(X,1);
X(:,num_zeros == 0) = [];

% division in training set and test set
X_training = X(1:length(X(:,1))*0.9,:);
Y_training = Y(1:length(X(:,1))*0.9,:);
X_test = X(length(X(:,1))*0.9:end,:);
Y_test = Y(length(X(:,1))*0.9:end,:);

% Matrix calculation
X_training = [ones(size(X_training,1),1) X_training];
X_test = [ones(size(X_test,1),1) X_test];
Y_tr = Y_training;

%% training

A = inv(transpose(X_training)*X_training + lambda *eye(size(X_training,2)))*transpose(X_training)*Y_training;

% Polynomial fitting (only one iteration)
Y_rs = X_training * A;

for ii = 1:size(Y_tr,2)
    [pol(ii,:),s,mupol(ii,:)] = polyfit(Y_rs(:,ii), Y_tr(:,ii),2);
end


%% Test

Y_res = X_test * A;
Y_fin = [];
time = linspace(1,(size(Y_res,1))*1/new_samp,size(Y_res,1));

% Polynomial application to the output of Wiener filter
for i=1:size(Y_res,2)
    Y_fin(:,i) = polyval(pol(i,:),Y_res(:,i),mupol(i,:));
end

% plot of the results
if f_plot
    figure  
    for k = 1:size(Y_test,2)
        sub(k) = subplot(size(Y_test,2),1,k);
        plot(time,Y_test(:,k))
        hold on
        plot(time,Y_res(:,k))
    end
end

%% evaluation
for k = 1:size(Y_test,2)
    R = corrcoef(Y_test(:,k),Y_fin(:,k));
    Test_perf(k) = R(1,2);
    MSE(k) = sum((Y_test(:,k) - Y_fin(:,k)).^2) / size(Y_test,1);
end

z = 1;
if find(num_zeros == 0,1,'first')
    while z < length(num_zeros)
        if num_zeros(z) == 0
            A = [A(1:z-1,:); repmat(0,1,size(A,2)); A(z:end,:)];
        end
        z = z+1;
    end
end


% explained VAF by the prediction
for k = 1:size(Y_test,2)
    VAF(k) = 1 - sum((Y_test(:,k) - Y_fin(:,k)).^2) / sum((Y_test(:,k) - mean(Y_test(:,k))).^2);
end


% plots for principal elements
% may have to be modified
if f_plot
    Wiener_plots;
end