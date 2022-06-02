#' Close-Ties matrix generator
#'
#' A renamed as.matrix function to take in an object that can be coerced into a matrix and returns it to the environment as a matrix ready for simulation.
#'@param close_matrix An object that can be coerced into a matrix. Values must be either 0 or 1 to denote connections.
#'@return Returns to the global environment a close adjacency matrix labeled adj_close.
#'@export

close_matrix = function(close_matrix){adj_close = as.matrix(close_matrix)
return(.GlobalEnv$adj_close <- adj_close)}

#' Weak-Ties matrix generator
#'
#' Generates an adjacency matrix for an environment where everyone is connected to everyone but themselves. Attempts to capture the social network at large.
#'@param ndata A dataset set up so that variable 1 is ID, variable 2 is behavior at time 0 (0-1 scale by fixed intervals), variable 3 is the close ties response, and variable 4 is the weak ties response. Total of 4 variables.
#'@return Returns to the global environment a weak-ties adjacency matrix labeled adj_weak.
#'@export

weak_matrix = function(ndata){ adj_weak = matrix(data = 1, nrow=nrow(ndata), ncol = nrow(ndata))
diag(adj_weak) = 0
return(.GlobalEnv$adj_weak <- adj_weak)}

#' Base simulation dataset generator
#'
#' @param dataset  A dataset set up so that variable 1 is ID, variable 2 is behavior at time 0 (0-1 scale by fixed intervals), variable 3 is the close ties response, and variable 4 is the weak ties response. Total of 4 variables.
#' @param close_matrix A unidirectional or bidirectional adjacency matrix. Values are either 0 or 1 to represent connections.
#' @return Returns a dataset, labeled ndata, to the global environment with simulation results.
#' @export

sim_dataset = function(dataset, close_matrix) {ndata = dataset
ndata$degree = NA
ndata$close_env = NA
ndata$weak_env = NA
ndata$time = 1
ndata$close_change = NA
ndata$weak_change = NA
ndata$joint_change = NA
ndata$bhv_temp = NA
names(ndata)  = c( "ID", "bhv_0", "close_ties", "weak_ties", "degree","close_env",
                   "weak_env", "time", "close_change", "weak_change", "joint_change", "bhv_temp" )
adj_close = close_matrix
ndata$degree = apply(adj_close,1,sum)
adj_weak = matrix(data = 1, nrow=nrow(ndata), ncol = nrow(ndata))
diag(adj_weak) = 0
ndata$degree = na_if(ndata$degree, 0)
ndata0 = ndata
return(.GlobalEnv$ndata <- ndata)}

#' Base simulation function
#'
#' @param dataset  A dataset set up so that variable 1 is ID, variable 2 is behavior at time 0 (0-1 scale by fixed intervals), variable 3 is the close ties response, and variable 4 is the weak ties response. Total of 4 variables.
#' @param close_matrix A unidirectional or bidirectional adjacency matrix. Values are either 0 or 1 to represent connections.
#' @param sim_max A value at which to limit the simulation. Warning, it may never reach that value if values do not improve from t to t+1.
#' @param bhv_increment The increments at which your behavior variable increases. It is a constant value.
#' @param CT_attribution_fraction The weight of close ties.
#' @param WT_attribution_fraction The weight of weak ties.
#' @return Returns a dataset, labeled ndata, to the global environment with simulation results.
#' @export

base_sim = function(dataset, close_matrix, sim_max, bhv_increment, CT_attribution_fraction, WT_attribution_fraction)
{
  ndata = dataset
  ndata$degree = NA
  ndata$close_env = NA
  ndata$weak_env = NA
  ndata$time = 1
  ndata$close_change = NA
  ndata$weak_change = NA
  ndata$joint_change = NA
  ndata$bhv_temp = NA
  names(ndata)  = c( "ID", "bhv_0", "close_ties", "weak_ties", "degree","close_env",
                     "weak_env", "time", "close_change", "weak_change", "joint_change", "bhv_temp" )
  adj_close = as.matrix(close_matrix)
  ndata$degree = apply(adj_close,1,sum)
  adj_weak = matrix(data = 1, nrow=nrow(ndata), ncol = nrow(ndata))
  diag(adj_weak) = 0
  ndata$degree = na_if(ndata$degree, 0)

  for (t in 1:sim_max) {
    ndata$bhv_temp = ndata[,paste0("bhv_",t-1)]
    ndata$close_env = as.vector((adj_close %*% ndata$bhv_temp) / ndata$degree)
    ndata$weak_env = as.vector((adj_weak %*% ndata$bhv_temp) / (nrow(ndata)-1))
    ndata$close_change = ifelse(ndata$close_ties <= ndata$close_env,1*CT_attribution_fraction,0)
    ndata$weak_change = ifelse(ndata$weak_ties <= ndata$weak_env,1*WT_attribution_fraction,0)
    ndata$joint_change = ndata$close_change + ndata$weak_change
    ndata[,paste0("bhv_",t)] = ifelse(ndata$bhv_temp + (bhv_increment)*ndata$joint_change >=1,
                                      1,ndata$bhv_temp + (bhv_increment)*ndata$joint_change)

    for (i in 1:length(ndata[,paste0("bhv_",t)])) {ifelse(is.na(ndata[i,paste0("bhv_",t)]),
                                                          ndata[i,paste0("bhv_",t)] <- ndata$bhv_0[i],
                                                          ndata[i,paste0("bhv_",t)] <- ndata[i,paste0("bhv_",t)]) }

    if ((sum(ndata[,paste0("bhv_",t)] ==ndata[,paste0("bhv_",t-1)]) == nrow(ndata)) == T)
      break #COMMENT: When no change, simulation is terminated
  }
  return(.GlobalEnv$ndata <- ndata)
}


#' Testing time simulation
#'
#' Without any public health interventions, how will behavior change with time?
#' This is a slight upgrade to base_sim that, in addition to producing 'ndata' dataset, also produces summary statistics.
#'
#' @param dataset  A dataset set up so that variable 1 is ID, variable 2 is behavior at time 0 (0-1 scale by fixed intervals), variable 3 is the close ties response, and variable 4 is the weak ties response. Total of 4 variables.
#' @param close_matrix A unidirectional or bidirectional adjacency matrix. Values are either 0 or 1 to represent connections.
#' @param sim_max A value at which to limit the simulation. Warning, it may never reach that value if values do not improve from t to t+1.
#' @param bhv_increment The increments at which your behavior variable increases. It is a constant value.
#' @param CT_attribution_fraction The weight of close ties.
#' @param WT_attribution_fraction The weight of weak ties.
#' @return Returns a dataset, labeled ndata, to the global environment with simulation results.
#' @export

tt= function(dataset, close_matrix, sim_max, bhv_increment,
                            CT_attribution_fraction, WT_attribution_fraction)
  {base_sim(dataset,
            close_matrix,
            sim_max,
            bhv_increment,
            CT_attribution_fraction,
            WT_attribution_fraction)
  cbind.data.frame(table(ndata[2]), table(ndata[ncol(ndata)]))
  return(summary(ndata[c(2,ncol(ndata))]))}

#' Riskiest people dataset
#'
#' Generates a dataset of individuals(and associated results) from testing time simulation.
#' These are the riskiest individuals at the end of tt simulation -- Time does not fix these individuals!
#'
#' @param tt_output Requires output dataset from tt function.
#' @param max_bhv Maximum behavior extraction value.
#' @return Returns a dataset, labeled riskiest_ppl, to the global environment with simulation results from tt.
#' @export

riskiest_ppl =  function(tt_output, max_bhv)
{ riskiest_ppl = tt_output[tt_output[length(tt_output)] <= max_bhv, ]
return(.GlobalEnv$riskiest_ppl <- riskiest_ppl)}

#' Permutation method
#'
#' Shuffles individuals around the network whilst maintaining structure and performs simulations.
#' Can be used to assess the distribution of simulation results.
#'
#' @param times The number of permutations to perform.
#' @param dataset  A dataset set up so that variable 1 is ID, variable 2 is behavior at time 0 (0-1 scale by fixed intervals), variable 3 is the close ties response, and variable 4 is the weak ties response. Total of 4 variables.
#' @param close_matrix A unidirectional or bidirectional adjacency matrix. Values are either 0 or 1 to represent connections.
#' @param sim_max A value at which to limit the simulation. Warning, it may never reach that value if values do not improve from t to t+1.
#' @param bhv_increment The increments at which your behavior variable increases. It is a constant value.
#' @param CT_attribution_fraction The weight of close ties.
#' @param WT_attribution_fraction The weight of weak ties.
#' @return Returns a dataset, labeled Perm_data, to the global environment with permutation results.
#' @export

permutate = function(times, dataset, close_matrix, sim_max, bhv_increment,
                     CT_attribution_fraction, WT_attribution_fraction )

{Perm_data= as.data.frame(1:times)
Perm_data$Results = NA
names(Perm_data) = c("Sim_Num", "Results")

for( m in 1:times){
  dataset[2:4] = dataset[sample(nrow(dataset)),2:4]
  base_sim(dataset, close_matrix, sim_max, bhv_increment,
           CT_attribution_fraction, WT_attribution_fraction)
  Perm_data[m, 2] = mean(ndata[, ncol(ndata)])
  print(paste0("Permutation #", m, " complete."))
}
return(.GlobalEnv$Perm_data <- Perm_data)}

#' Individual approach simulation
#'
#' Modifies, in order, the behavior of an individual by one behavior increment and runs the simulation thereafter. It then stores the results, resets, moves on to the next person, and continues until every individual's impact has been recorded.
#'
#' @param dataset  A dataset set up so that variable 1 is ID, variable 2 is behavior at time 0 (0-1 scale by fixed intervals), variable 3 is the close ties response, and variable 4 is the weak ties response. Total of 4 variables.
#' @param close_matrix A unidirectional or bidirectional adjacency matrix. Values are either 0 or 1 to represent connections.
#' @param sim_max A value at which to limit the simulation. Warning, it may never reach that value if values do not improve from t to t+1.
#' @param bhv_increment The increments at which your behavior variable increases. It is a constant value.
#' @param CT_attribution_fraction The weight of close ties.
#' @param WT_attribution_fraction The weight of weak ties.
#' @return Returns a dataset, labeled Ind_data, to the global environment with simulation results.
#' @export

ind_approach = function(dataset, close_matrix, sim_max, bhv_increment,
                        CT_attribution_fraction, WT_attribution_fraction )

{x = nrow(dataset)
Ind_data= as.data.frame(1:x)
Ind_data$Results = NA
names(Ind_data) = c("Sim_ID_Changed", "Results")
datasettemp = dataset
for( m in 1:x){
  datasettemp = dataset
  ifelse( datasettemp[m, 2] < 1, datasettemp[m, 2] <- datasettemp[m, 2]+ bhv_increment, datasettemp[m, 2] <- datasettemp[m, 2])
  base_sim(datasettemp, close_matrix, sim_max, bhv_increment,
           CT_attribution_fraction, WT_attribution_fraction)
  Ind_data[m, 2] = mean(ndata[, ncol(ndata)])
  print(paste0("Outcome of individual ", m, " behavior recorded."))
}
return(.GlobalEnv$Ind_data <- Ind_data)}

#' Most influential individuals
#'
#' Generates a dataset of the most influential individuals from the individual approach(ind_approach function).
#'
#' @param ind_data_output  Requires output dataset from ind_approach function.
#' @param testing_time_avg The average value from the final time period of the tt function. Used to standardize results.
#' @return Returns a dataset, labeled influential_inds, to the global environment with simulation results.
#' @export

Order_of_influence = function(ind_data_output, testing_time_avg)

{ind_data_output$Change = ind_data_output$Results - testing_time_avg
influential_inds <- ind_data_output[order(-ind_data_output$Change),]
row.names(influential_inds) = NULL
return(.GlobalEnv$influential_inds <- influential_inds)}

#' High risk approach
#'
#' Simulate behavior change among high risk individuals (riskiest_ppl output).
#'
#' @param num_of_ind The number of individuals to randomly select from riskiest people dataset.
#' @param coercion_factor Behavior increment to boost behavior by. Regardless of number, individuals will be capped at 1.
#' @param riskiest_ppl_dataset Requires output dataset from riskiest_ppl function.
#' @param dataset  A dataset set up so that variable 1 is ID, variable 2 is behavior at time 0 (0-1 scale by fixed intervals), variable 3 is the close ties response, and variable 4 is the weak ties response. Total of 4 variables.
#' @param close_matrix A unidirectional or bidirectional adjacency matrix. Values are either 0 or 1 to represent connections.
#' @param sim_max A value at which to limit the simulation. Warning, it may never reach that value if values do not improve from t to t+1.
#' @param bhv_increment The increments at which your behavior variable increases. It is a constant value.
#' @param CT_attribution_fraction The weight of close ties.
#' @param WT_attribution_fraction The weight of weak ties.
#' @return Returns a dataset, labeled ndata, to the global environment with simulation results.
#' @export

high_risk_approach = function(num_of_ind, coercion_factor, riskiest_ppl_dataset, dataset, close_matrix, sim_max, bhv_increment,
                              CT_attribution_fraction, WT_attribution_fraction)
{ data_temp = dataset
x = num_of_ind #Number of individuals to select
y = coercion_factor #Behavior coercion factor
outcome = sort(sample(riskiest_ppl_dataset$ID, x, replace = FALSE))
print(paste0("The following indviduals' behavior has been coerced by a factor of ", coercion_factor, "."))
print(outcome)
data_temp[data_temp$ID %in% outcome,2] = data_temp[data_temp$ID %in% outcome,2] + y
data_temp[data_temp[2] >= 1, 2] = 1
base_sim(data_temp, close_matrix, sim_max, bhv_increment,
         CT_attribution_fraction, WT_attribution_fraction)
}

#' Population approach
#'
#' Simulate behavior change among randomly selected individuals. Excludes those already at 1.
#'
#' @param num_of_ind The number of individuals to randomly select from the population.
#' @param coercion_factor Behavior increment to boost behavior by. Regardless of number, individuals will be capped at 1.
#' @param dataset  A dataset set up so that variable 1 is ID, variable 2 is behavior at time 0 (0-1 scale by fixed intervals), variable 3 is the close ties response, and variable 4 is the weak ties response. Total of 4 variables.
#' @param close_matrix A unidirectional or bidirectional adjacency matrix. Values are either 0 or 1 to represent connections.
#' @param sim_max A value at which to limit the simulation. Warning, it may never reach that value if values do not improve from t to t+1.
#' @param bhv_increment The increments at which your behavior variable increases. It is a constant value.
#' @param CT_attribution_fraction The weight of close ties.
#' @param WT_attribution_fraction The weight of weak ties.
#' @return Returns a dataset, labeled ndata, to the global environment with simulation results.
#' @export

population_approach = function(num_of_ind, coercion_factor, dataset, close_matrix, sim_max, bhv_increment,
                               CT_attribution_fraction, WT_attribution_fraction)
{ data_temp = dataset
x = num_of_ind #Number of individuals to select
y = coercion_factor #Behavior coercion factor

out = data_temp[!data_temp[2] >= 1,"ID" ]
randos = sort(sample(out, x, replace = FALSE))

print(paste0("The following indviduals' behavior has been coerced by a factor of ", coercion_factor, "."))
print(randos)

data_temp[data_temp[,1] %in% randos, 2] = data_temp[data_temp[,1] %in% randos, 2] + y
data_temp[data_temp[2] >= 1, 2] = 1

base_sim(data_temp, close_matrix, sim_max, bhv_increment,
         CT_attribution_fraction, WT_attribution_fraction)
}

#' Influential approach
#'
#' Simulate behavior change among the most influential individuals. This is determined by the ranked output from Order_of_influence function.
#'
#' @param num_of_ind The number of individuals to select from the list of most influential individuals. Selects individuals in order. ex. 5 would be the top most influential individuals to coerce.
#' @param coercion_factor Behavior increment to boost behavior by. Regardless of number, individuals will be capped at 1.
#' @param ordered_influence_dataset  Dataset from Order_of_influence function; labeled influential_inds.
#' @param dataset  A dataset set up so that variable 1 is ID, variable 2 is behavior at time 0 (0-1 scale by fixed intervals), variable 3 is the close ties response, and variable 4 is the weak ties response. Total of 4 variables.
#' @param close_matrix A unidirectional or bidirectional adjacency matrix. Values are either 0 or 1 to represent connections.
#' @param sim_max A value at which to limit the simulation. Warning, it may never reach that value if values do not improve from t to t+1.
#' @param bhv_increment The increments at which your behavior variable increases. It is a constant value.
#' @param CT_attribution_fraction The weight of close ties.
#' @param WT_attribution_fraction The weight of weak ties.
#' @return Returns a dataset, labeled ndata, to the global environment with simulation results.
#' @export

network_approach = function(num_of_ind, coercion_factor, ordered_influence_dataset, dataset, close_matrix, sim_max, bhv_increment,
                            CT_attribution_fraction, WT_attribution_fraction)
{ data_temp = dataset

x = num_of_ind #Number of individuals to select
y = coercion_factor #Behavior coercion factor

out = ordered_influence_dataset[ 1:x,]

print(paste0("The following indviduals' behavior has been coerced by a factor of ", coercion_factor, "."))
print(as.vector(out[,1]))

data_temp[data_temp[,1] %in% out[,1], 2] = data_temp[data_temp[,1] %in% out[,1], 2] + y
data_temp[data_temp[2] >= 1, 2] = 1

base_sim(data_temp, close_matrix, sim_max, bhv_increment,
         CT_attribution_fraction, WT_attribution_fraction)
}

#' Rough cost benefit analysis
#'
#' Produces a rough cost benefit analysis. Utilizes the influential individuals dataset as they are the most likely to influence other.
#'
#' @param bin_size The number of units to be spent. Simulation stops if maximum results are obtained prior to utilizing all units. One unit is equal to one behavior increment.
#' @param ordered_influence_dataset  Dataset from Order_of_influence function; labeled influential_inds.
#' @param dataset  A dataset set up so that variable 1 is ID, variable 2 is behavior at time 0 (0-1 scale by fixed intervals), variable 3 is the close ties response, and variable 4 is the weak ties response. Total of 4 variables.
#' @param close_matrix A unidirectional or bidirectional adjacency matrix. Values are either 0 or 1 to represent connections.
#' @param sim_max A value at which to limit the simulation. Warning, it may never reach that value if values do not improve from t to t+1.
#' @param bhv_increment The increments at which your behavior variable increases. It is a constant value.
#' @param CT_attribution_fraction The weight of close ties.
#' @param WT_attribution_fraction The weight of weak ties.
#' @return Returns a dataset, labeled bin_data, to the global environment with results.
#' @export



rough_cba = function(bin_size, ordered_influence_dataset,
                     dataset, close_matrix, sim_max, bhv_increment,
                     CT_attribution_fraction, WT_attribution_fraction)

{initial_bin_size = 1:bin_size
final_bin_size = NA
individuals_reformed = NA
result = NA
bin_data = cbind.data.frame(initial_bin_size, final_bin_size, individuals_reformed, result)

#run it

for (n in 1:bin_size){
  bin_data$initial_bin_size[n] = n

  temp_dataset = dataset

  bin = n

  for (k in 1:173){
    person = ordered_influence_dataset[ k,1]
    bhv = temp_dataset[temp_dataset[,1] == person, 2]
    used = (1 - bhv)/(bhv_increment)
    temp_dataset[temp_dataset[,1] == person, 2] <- bhv + used*(bhv_increment)
    temp_bin = bin - used
    bin = temp_bin

    # print(paste("Individual", person, "Reformed."))

    if ((bin) < (1 - temp_dataset[temp_dataset[,1] == ordered_influence_dataset[ifelse(k < 173,k+1,k),1], 2])/(1/3))
    {print(paste("Bin value is", round(bin,5), ".", k, "total individuals reformed.",
                 "Coersion of individual",
                 ordered_influence_dataset[ k+1,"Sim_ID_Changed"], "requires additional units. Please reload."))
      break}
  }
  bin_data$final_bin_size[n] = bin


  bin_data$individuals_reformed[n] = k

  base_sim(temp_dataset,
           close_matrix,
           sim_max,
           bhv_increment,
           CT_attribution_fraction,
           WT_attribution_fraction)

  bin_data$result[n] = mean(ndata[,ncol(ndata)])

  if (bin_data$result[n] == 1) {
    print(paste("Bin value is maximized at", n))

    bin_data = bin_data[complete.cases(bin_data),]

    return(.GlobalEnv$bin_data <- bin_data)
    break}
}}

