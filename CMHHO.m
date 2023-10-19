% Developed in MATLAB R2013b
% Harris's hawk optimizer
% Dimension decided Harris hawks optimization with Gaussian mutation: Balance analysis and diversity patterns
% https://doi.org/10.1016/j.knosys.2020.106425

% T: maximum iterations, N: populatoin size, CNVG: Convergence curve
%N 为30，D为30； MaxFES为最大评估次数设为300000；每次实验重复30次
function [Rabbit_Location, CNVG]=RHHO_fes(N,MaxFES,lb,ub,dim,fobj)

disp('HHO is now tackling your problem')
tic
% initialize the location and Energy of the rabbit
Rabbit_Location=zeros(1,dim);
Rabbit_Energy=inf;
% N行1列全为正无穷的矩阵

%Initialize the locations of Harris' hawks
X=initialization(N,dim,ub,lb);
d=dim;
CNVG=[];
fes=0; % FES counter

t=0; % Loop counter

while fes<MaxFES
    for i=1:size(X,1)
        % Check boundries
        FU=X(i,:)>ub;FL=X(i,:)<lb;X(i,:)=(X(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
        % fitness of locations
        fitness(i)=fobj(X(i,:));
        fes=fes+1;
        % Update the location of Rabbit
        if fitness(i)<Rabbit_Energy
            Rabbit_Energy=fitness(i);
            Rabbit_Location=X(i,:);
        end
    end
    
    c=2*(1-(fes/MaxFES)); % factor to show the decreaing energy of rabbit
    % Update the location of Harris' hawks
    for i=1:size(X,1)
        Escaping_Energy=c*(2*rand()-1);  % escaping energy of rabbit
        
        if abs(Escaping_Energy)>=1
            %% Exploration:
            % Harris' hawks search for the rabbit, observing the area to find a rabbit
            % Harris' hawks perch randomly based on 2 strategy:
            
            r=rand();
            rand_Hawk_index = floor(N*rand()+1);
            X_rand = X(rand_Hawk_index, :);
            if r<0.5
                % perch based on other family members
                X(i,:)=X_rand-rand()*abs(X_rand-2*rand()*X(i,:));
            elseif r>0.5
                % perch on a random tall tree (random site inside group's home range)
                X(i,:)=(Rabbit_Location(1,:)-mean(X))-rand()*((ub-lb)*rand+lb);
            end
            
        elseif abs(Escaping_Energy)<1
            %% Exploitation:
            % Attacking the rabbit using 4 strategies regarding the behavior of the rabbit
            
            %% phase 1: surprise pounce (seven kills): hawks go to kill when rabbit is surprised
            % surprise pounce (seven kills): multiple, short rapid dives by different hawks
            
            r=rand(); % probablity of each event % r>0.5 : rabbit cannot escape -- r<0.5 : rabbit can escape
            
            if r>0.5 && abs(Escaping_Energy)<0.5 % Hard besiege % When rabbit cannot or do not try to escape
                X(i,:)=(Rabbit_Location)-Escaping_Energy*abs(Rabbit_Location-X(i,:));
            end
            
            if r>0.5 && abs(Escaping_Energy)>0.5  % Soft besiege % When rabbit cannot or do not try to escape
                Jump_strength=2*(1-rand()); % random jump strength of the rabbit
                X(i,:)=(Rabbit_Location-X(i,:))-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:));
            end
            
            %% phase 2: performing team rapid dives (leapfrog movements) when prey escapes
            if r<0.5 && abs(Escaping_Energy)>0.5 % Soft besiege % rabbit try to escape by many zigzag deceptive motions
                
                Jump_strength=2*(1-rand());
                X1=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:));
                
                fes=fes+1;
                if fobj(X1)<fitness(i) % improved move?
                    X(i,:)=X1;
                else % hawks perform levy-based short rapid dives around the rabbit
                    X2=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:))+rand(1,dim).*Levy(dim);
                    fes=fes+1;
                    if fobj(X2)<fitness(i) % improved move?
                        X(i,:)=X2;
                    end
                end
            end
            
            if r<0.5 && abs(Escaping_Energy)<0.5 % Hard besiege % rabbit try to escape by many zigzag deceptive motions
                % hawks try to decrease their average location with the rabbit
                Jump_strength=2*(1-rand());
                X1=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-mean(X));
                
                fes=fes+1;
                if fobj(X1)<fitness(i) % improved move?
                    X(i,:)=X1;
                else % Perform levy-based short rapid dives around the rabbit
                    X2=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-mean(X))+rand(1,dim).*Levy(dim);
                    fes=fes+1;
                    if fobj(X2)<fitness(i) % improved move?
                        X(i,:)=X2;
                    end
                end
            end
           
        end
        
           %% 柯西变异更新位置
           A = randi([1,N]);  %随机选择N个种群随机一个
           B = randi([1,N]);
           %%变异 mutation
           newX = Rabbit_Location-cauchy(0,0.5,d).*(X(A,:)-X(B,:)); %柯西变异和随机个体更新位置
           fes=fes+1;
           if fobj(newX)<fitness(i)
            X(i,:)=newX;
           end
             %% 
    end
   
    t=t+1;
    CNVG(t)=Rabbit_Energy;
    %Print the progress every 50 iterations
%     if mod(t,50)==0
%         display(['At iteration ', num2str(t), ' the best fitness is ', num2str(Rabbit_Energy)]);
%     end
end
toc
end

function o=Levy(d)
beta=1.5;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;v=randn(1,d);step=u./abs(v).^(1/beta);
o=step;
end