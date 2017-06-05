function meanSim = bootstrapSim (sim1,sim2,itr)

sim1 = sim1-mean(sim1);
sim2 = sim2-mean(sim2);
n1 = length(sim1);
n2 = length(sim2);
meanSim = zeros(itr,2);
for j = 1:itr
     p1 = randsample(n1,n1,true);
     p2 = randsample(n2,n2,true);
     boot_sim1 = sim1(p1);
     boot_sim2 = sim2(p2);
     meanSim(j,:) = [mean(boot_sim1) mean(boot_sim2)];
end