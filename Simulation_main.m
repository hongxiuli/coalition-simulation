%main function for the simulation
%parameter values are imported from the "parameter.xls"
betapos = -beta;
phi = vul;
n = length(alpha);

%index is a matrix containing all possible combinations of the elements of vector [1:10] taken s at a time.
for s = 1:n
  index = nchoosek(1:n,s); % choose s from 1:n sequence
  stable_coal(s).size = s;
  [intcoal,extcoal,coalition,E_S,E,W] = IEA_adapt_func(alpha,betapos,phi, index, n);
  stable_coal(s).coalition = coalition;
  stable_coal(s).intcoal = intcoal;
  stable_coal(s).extcoal = extcoal;
  stable_coal(s).ems_IEA = E_S;
  stable_coal(s).ems = E;
  stable_coal(s).welfare = W;
end


coal = [2,4,5,6];
%coal =[10];
O = ones(10,1);
O(coal) = 0;
noncoal = find(O); %keep all non-zero elements in O
r_graph = [0.5:0.01:1]; %row vector
r_S = 1;
gamma_graph = zeros(length(noncoal),length(r_graph));
gamma_temp_i = zeros(length(noncoal),length(r_graph));
for m = 1:length(r_graph)
  r = r_graph(m);
  [gamma_graph(:,m),gamma_temp_i(:,m)] = incentive_adapt(alpha,betapos,phi,noncoal,coal,r,r_S);
end

legendCell = cellstr(num2str(noncoal, 'i=%-d'));
%axis([0,1,-1.2*10^10,1.8*10^8])   
figure
plot(r_graph,gamma_graph(1,:),':', r_graph,gamma_graph(2,:),'-.',r_graph,gamma_graph(3,:),'-o',r_graph,gamma_graph(4,:),'--',r_graph,gamma_graph(5,:),'-', r_graph,gamma_graph(6,:), '-x',...
    'LineWidth',1.5,...
    'MarkerSize',3)
xlabel('Adoption cost r') % x-axis label
ylabel('Free-riding incentive of non-members') % y-axis label
legend(legendCell,'Location','southeast')


