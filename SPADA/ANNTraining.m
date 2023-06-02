function output = ANNTraining(filename, filepath)
% Data Processing and train an ANN
% Remove first principal strain outliners 
% output: angular deflection

DataSetFile = fullfile(filepath, filename);
load(DataSetFile);
N_sample = length(sample_no);
S = size(result_value);
P_step = S(2);
X = [];
y = [];
Outliner = 0;

for i = 1: N_sample
    Outliner_bool = 0;
    x = zeros(5,P_step);
    for j = 1:P_step
        if result_value(3,j,i) > 3 % critial strain value
            Outliner = Outliner + 1;
            Outliner_bool = 1;
            break
        end
            x(1:4,j) = variable_value(:,1,i);
    end
    if Outliner_bool == 1
        continue
    else
        x(5,:) = result_value(1,:,i);
        X = [X x];
        r = result_value(2,:,i);
        y = [y r];
    end
end

% ANN training
net = feedforwardnet;
net=configure(net,X,y);
trainFcn = 'trainbr';
net = feedforwardnet([20 20],trainFcn);
net.divideFcn = 'dividerand';
net.divideParam.trainRatio = 80/100;
net.divideParam.valRatio = 0/100;
net.divideParam.testRatio = 20/100;
net.performFcn = 'mse';
net.performParam.normalization = 'standard';
net.trainParam.epochs = 1000;

% Set the output function to update the window position
net.trainParam.OutputFcn = @(info) plotWindowPosition(info, windowPosition);
net = train(net, X, y);

% save net
outputfile = fullfile(pwd, "NewANN.mat");
save(outputfile, 'net');
output = "NewANN.mat";
end

function plotWindowPosition(info, windowPosition)
    % Retrieve the handle of the training progress plot
    hPlot = gcf;
    % Get screen size
    screenSize = get(groot, 'ScreenSize');
    screenWidth = screenSize(3);
    screenHeight = screenSize(4);
    figWidth = 800;
    figHeight = 600;
  
    % Set the window position using the specified coordinates
    set(hPlot, 'Position', [(screenWidth - figWidth)/2, (screenHeight - figHeight)/2, figWidth, figHeight]);  % Adjust width and height as desired
end
