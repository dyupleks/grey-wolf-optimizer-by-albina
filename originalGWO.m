clc;
clear;
close all;

%%
tic
hf = figure; %Open figure and keep handle
hf=colordef(hf,'white'); %Set color scheme
hf.Color='w'; %Set background color of figure window

%%
CostFunction = @(x) Sphere(x);
dim = 2;
varSize = [1 dim];

ub = 5;  %upper bound
lb = 0;  %lower bound

nPop = 10;
maxIt = 20;

empty_wolf.Position = [];
empty_wolf.Cost = [];
wolfs = repmat(empty_wolf, nPop, 1);



for i=1:nPop
    wolfs(i).Position = unifrnd(lb, ub, varSize);
    wolfs(i).Cost = CostFunction(wolfs(i).Position);
    wolfs(i).Best.Cost = wolfs(i).Cost;
    wolfs(i).Best.Position = wolfs(i).Position; 
end

alpha.Cost = inf;
beta.Cost = inf;
gamma.Cost = inf;


% run only once in initialization
for i=1:nPop
   if wolfs(i).Best.Cost < alpha.Cost
      gamma = beta;
      beta = alpha;
      alpha = wolfs(i).Best;
   
   elseif wolfs(i).Best.Cost < beta.Cost
       gamma = beta;
       beta = wolfs(i).Best;
       
   elseif wolfs(i).Best.Cost < gamma.Cost
       gamma = wolfs(i).Best;
   end    
end

a = 2;
adamp = a/maxIt;
l=1;
%array to hold the alpha costs during whole iterations
alphaCosts = zeros(maxIt,1);

map = robotics.BinaryOccupancyGrid (5, 5);
show(map);
hold on
axis
xlabel('x');
ylabel('y');
% xlim([-4 6])
% ylim([-4 6])
hold on

%% plotting parameters

width = 5;     % Width in inches
height = 5;    % Height in inches
alw = 0.95;    % AxesLineWidth
fsz = 17;      % Fontsize
lw = 3;      % LineWidth
msz = 8;       % MarkerSize

% pos = get(gcf, 'Position');
% set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
%%


%Give colors
grid off

colors = lines(nPop);
f(1 : nPop) = animatedline();
h(1 : nPop) = animatedline();

for i=1:nPop
   f(i) = animatedline('Color', colors(i, :),'LineWidth',lw);
   h(i) = animatedline('Marker','o','MarkerFaceColor',colors(i, :),'MarkerSize',12);
end

%%
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
%%


%Main loop

for it=1:maxIt
    for i=1:nPop 
        x(i, it) = wolfs(i).Best.Position(1,1);
        y(i, it) = wolfs(i).Best.Position(1,2);
        
        
    end 
   
    for j=1:100
    
        for i=1:nPop 
            if it == 1 
                xlsi = x(i, it);
                addpoints(h(i),xlsi,y(i, it));
                addpoints(f(i),xlsi,y(i, it));   
                
                
                hold on            
            else               
                xlsi = linspace(x(i, it-1),x(i, it));
                ylsi = linspace(y(i, it-1),y(i, it));
                addpoints(h(i),xlsi(j),ylsi(j));
                addpoints(f(i),xlsi(j),ylsi(j));
                
                 drawnow
%                 drawnow limitrate   
                clearpoints(h(i))
                addpoints(h(i),xlsi(j),ylsi(j));   

            end
        end
    end
    xlsi = x(i, it);
    
    addpoints(h(i),xlsi,y(i, it));

    axis 'auto xy';
    title(['Iteration = ', num2str(it)]);

    A1 = 2*a*rand(varSize) - a;
    C1 = 2*rand(varSize);
    A2 = 2*a*rand(varSize) - a;
    C2 = 2*rand(varSize);
    A3 = 2*a*rand(varSize) - a;
    C3 = 2*rand(varSize);
    
    for i=1:nPop
        Dalpha = abs(C1.*alpha.Position - wolfs(i).Position);
        X1 = alpha.Position - A1.*Dalpha;
        Dbeta = abs(C2.*beta.Position - wolfs(i).Position);
        X2 = beta.Position - A2.*Dbeta;
        Dgamma = abs(C3.*gamma.Position - wolfs(i).Position);
        X3 = gamma.Position - A3.*Dgamma;
        
        wolfs(i).Position = (X1 + X2 + X3)/3;
        wolfs(i).Cost = CostFunction(wolfs(i).Position);
        wolfs(i).Best.Cost = wolfs(i).Cost;
        wolfs(i).Best.Position = wolfs(i).Position;         
    end
    
    alphaCosts(it) = alpha.Cost;
    for i=1:nPop
        if wolfs(i).Best.Cost < alpha.Cost
            gamma = beta;
            beta = alpha;
            alpha = wolfs(i).Best;
   
        elseif wolfs(i).Best.Cost < beta.Cost
            gamma = beta;
            beta = wolfs(i).Best;
       
        elseif wolfs(i).Best.Cost < gamma.Cost
            gamma = wolfs(i).Best;
        end    
    end 
    if it==1
       a;
    else
        a = a - adamp;
    end
   disp(['Iteration ' num2str(it) ': Alpha Cost = ' num2str(alphaCosts(it))]);
   disp([': a = ' num2str(a)]);
   

end

%% below code lines are only for visualization (no GWO)
width = 5;     % Width in inches
height = 5;
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
% 
figure(2);

hf2 = figure(2); %Open figure and keep handle
hf2=colordef(hf2,'white'); %Set color scheme
hf2.Color='w'; %Set background color of figure window
xlim([0,maxIt])
set(gca,'XTick',alphaCosts(maxIt):alphaCosts(1));
set(gca,'YTick',0:maxIt);
semilogy(alphaCosts, 'LineWidth', 4);
xlabel('Iterations'); 
ylabel('Alpha cost ');

set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

title('Alpha wolf optimization');


%     % Save the file as PNG
print(hf2,'alphacost','-dpng','-r200');

toc
