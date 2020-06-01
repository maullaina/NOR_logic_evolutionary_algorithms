% Continuous Genetic Algorithm 

%Aina maull y Marian Iglesias

% I Setup the GA___________________________________________________________
clear all

% objective function
ff= 'nueva_moap_0706'; 

% variable limits, rangos, el gen no puede tener cualquier valor
varhi= 15;  
varlo= 1; 

% II Stopping criteria_____________________________________________________
maxit=400;  % max number of iterations
mincost= -9999999; % minimum cost

% III GA parameters________________________________________________________
popsize =2500; % set population size
selection = 0.7;% fraction of population kept
mutrate_1 = 0.15;
mutrate_2 = 0.4; % set mutation rate
keep=floor(selection*popsize); % #population memberst that survive
M=ceil((popsize-keep)/2); % number of matings
% %pesos
 w = ones(1,8); %inicialmente todos los pesos valen 1
%  y = zeros(1,popsize);
%  ev = zeros(1,popsize);
 
 

% Create the initial population____________________________________________
iga=0; % generation counter initialized
cost = zeros(1,popsize);

C = {};%creem un Cell Array de dimensions desconegudes que anirem omplint
for a = 1:popsize %iterem fins completar tota la població
    O =randi([varlo,varhi]);%definim un valor random per les columnes
    Mi=round(rand(O+2,O));%matriu base per cada individu
    for x = 1:(O+2)%li sumem 2 perque hi ha 2 files mes
        for y = 1:O
            if  x > (y+2)%fem que no hi hagi zeros ala diagonal de baix
                 Mi(x,y) = 0;
            end   
        end
    end
    C{a}=Mi;
end

%Funcion FITNESS___________________________________________________________

for a=1:popsize
   [cost(a),w]=feval(ff,C{a},w);%feval es una instruccion para ejecutar una funcion con los parametros seguidos
end
[cost,p] =sort(cost);% min cost in element 1

%seleccio__________________________________________________________________
Sel = {};
Sel = C;
for i=1:popsize
    C{i} = Sel{p(i)};
end 
minc(1)=min(cost); % minc contains the minimum of population
meanc(1)=mean(cost); % meanc contains average fitness of population

while iga < maxit
    iga=iga+1 % increments generation counter
    %mutamos
    for i=1:keep
        CFit{i} = C{i};
    end 
R= round((6/10)*(popsize-keep));%fraccio que es recombinin (triar un num en concret)
Rec={};
Rec_Mut={};
Mut={};
for i = 1: R %tots els que es recombinen muten -- podria millorarse pero he de pensa com treure matrius nules dels arrays
    num1=randi([1,keep]);
    num2=randi([1,popsize]);
    s1= C{num1};%seleccio de dos fills random en concret que es recombinaran
    s2= C{num2};
    s1z= size(s1);
    s2z=size(s2);
      if s1z(1)>s2z(1)
        s_nw=zeros(s2z(1),s2z(2));
        for a = 1 : s2z(1)
            for b=1:s2z(2)
                s_nw(a,b)=s1(a,b);
            end
        end
        s1_nw=s_nw;
        s2_nw=s2;
      
      else
        s_nw=zeros(s1z(1),s1z(2));
        for a = 1 : s1z(1)
            for b=1:s1z(2)
                s_nw(a,b)=s2(a,b);  
            end
        end
        s2_nw=s_nw;
        s1_nw=s1;
      end
      
    sz_nw=size(s1_nw);
    rd=randi([1,sz_nw(2)]);
    s = [s1_nw(:,1:rd),s2_nw(:,rd+1:sz_nw(2))];
    szs=size(s);
    Rec{i}=s;
    m = rand(szs(1),szs(2)) < mutrate_1;
       for x = 1: szs(1);
          for y = 1:szs(2);
                  if m(x,y) == 1 && x >(y+2) %recorro la matriz de mutacion y miro que sea circuito valido
                         m(x,y) = 0;
                        end
             end
       end
         Rec_Mut{i} = mod(Rec{i} + m,2); %nueva poblacion mutada 

end
%petita porcio de mutacio -- anivellada amb la de recombinacio
T=round(0.8*(popsize-(keep+R)));
for j=1:T
num1=randi([1,5]);
s1= C{num1};
s1z=size(s1);
   if rand(1) < mutrate_2
        m = rand(s1z(1),s1z(2)) < mutrate_1;
       for x = 1: s1z(1)
          for y = 1:s1z(2)
                  if m(x,y) == (1 && x >(y+2)) %recorro la matriz de mutacion y miro que sea circuito valido
                         m(x,y) = 0;
                  end
          end
       end
         Mut{j} = mod(s1+ m,2); %nueva poblacion mutada
   else
            Ccf = C{num1};
            c = randi(1,s1z(2));
            d = c+2;
            Ccf = [Ccf(:,1:c),Ccf(:,c+2:s1z(2))]; %eliminamos la columna
            Ccf = [Ccf(1:d,:);Ccf(d+2:end,:)];
            Mut{j} = Ccf;
   end
end

%inmigracio
I= popsize - (keep + R + T);
Inm={};
for a = 1:I%recorro tota la població
    O =randi([varlo,varhi]);
    Mi=round(rand(O+2,O));%matriu base per cada individu
    for x = 1:(O+2)
        for y = 1:O
            if  x > (y+2)
                 Mi(x,y) = 0;
            end   
        end
    end
    Inm{a}=Mi;
end

C=horzcat(CFit, Rec_Mut, Mut, Inm);

    % The new offspring and mutated chromosomes are evaluated
    % ---------------------------------------------------
    for a=1:popsize
        [cost(a),w]=feval(ff,C{a},w);
%        [x(a),y(a),cost(a)]=feval(ff,C{a},x(a),y(a),cost(a));;
    end
    % Sort the costs and associated parameters
    % ---------------------------------------------------
    [cost, p] =sort(cost);
    Sel = {};
    Sel = C;
    for i=1:popsize
        C{i} = Sel{p(i)};
    end 
    best = C{1}
    mid = C{round(popsize/2)}
    % Do statistics for a single nonaveraging run
    % ---------------------------------------------------
    minc(iga+1)=min(cost);
    meanc(iga+1)=mean(cost);
    
    % Stopping criteria
    if iga>maxit || cost(1)<mincost
        break 
    end
    [iga cost(1)]; 
end
day=clock; 
disp(datestr(datenum(day(1),day(2),day(3),day(4),day(5), day(6)),0))
disp(['optimized function is ' ff])
format short g
disp(['popsize = ' num2str(popsize) ' mutrate_1 = ' num2str(mutrate_1) ' # mutrate_2 = ' num2str(mutrate_2)]) 
disp(['#generations=' num2str(iga) ' best cost=' num2str(cost(1))])
disp('best solution') 
disp(C{1}) 
disp('continuous genetic algorithm')
figure(1)
iters=0:length(minc)-1;
plot(iters,minc,iters,meanc,'-');
xlabel('generation'); ylabel('fitness');
% surface(x,y,ev);

     
 