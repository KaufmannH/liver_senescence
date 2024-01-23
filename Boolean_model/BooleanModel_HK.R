' This script can generate a Boolean Model that shows the interactions between 
senescent cells, and its microenvironment that consists of T lymphocytes (T_Lymph),
pro-inflammatory macrophages (Macro), anti-inflammatory macrophages (Supr_macro), 
resistant senescent cells   .'


# import packages
library(BoolNet)
library(here)
library(gplots)

# import rules and and load model
rules <- here("bool_v15.txt")
suppressWarnings(model <- loadNetwork(rules))

# rename genes to nodes
model$nodes <- model$genes

# option to visualize the connections of the model
# graph <- plotNetworkWiring(model)


## SETTINGS ##

# set initial states of nodes (0 or 1)
init.state <- generateState(model, specs = c(
  "anti_inf_sasp" = 0,
  "anti_inf_macro" = 0,
  "inf_macro" = 0,
  "res_sc" = 0,
  "sc" = 1,
  "T_cell" = 0
))

# set delays
# fill in the number of iterations after which the node can switch on
delays <- c("anti_inf_sasp" = 10,  "anti_inf_macro" = 10, "inf_macro" = 5,
            "res_sc" = 30, "sc" = 0, "T_cell" = 15)

# set off switches
# set the negative number of iterations after which the node should be switched off
off.switch <- c("anti_inf_sasp" = -1000, "anti_inf_macro" = -1000, "inf_macro" = -1000,
                "res_sc" = -1000, "sc" = -1000, "T_cell" = -1000)



## FUNCTIONS ##

Interventions <- function(trans, timesteps.counter, off, remaining.delays, model, init.state) {
'This function takes the current timestep, checks the requirements and inhibits certain nodes
  and returns the changed model.'
  
  # copy of the model and init.state to track knocked-out genes
  updated.model <- model 
  updated.state <- init.state
  
  for (g in 1:length(updated.model$nodes)) {
    off[g] <- off[g] + 1
    
    # calculate delays
    if (remaining.delays[g] > 0) { 
      gene.name <- names(remaining.delays[g])
      gene.name.vector <- c(gene.name) 
      updated.model <- fixGenes(updated.model, gene.name.vector, 0)
      updated.state[gene.name] <- 0 
      init.state <- updated.state
      remaining.delays[g] <- remaining.delays[g] - 1 
    }
    
    # Switch gene off
    if (off[g] > 0) { 
      gene.name <- names(off[g])
      gene.name.vector <- c(gene.name) 
      updated.model <- fixGenes(updated.model, gene.name.vector, 0)
      updated.state[gene.name] <- 0
      init.state <- updated.state
    }
    
    # Resurrection stop
    # stop senescent cell from switching on again after it has been off once.
    if (updated.model$nodes[g] == "sc"  && trans[g] == 0 && timesteps.counter > 1){
      sc.name.vector <- c("sc") 
      updated.model <- fixGenes(updated.model, sc.name.vector, 0)
      updated.state["sc"] <- 0
      init.state <- updated.state
    }
  }
  return(list(off = off, updated.model = updated.model, init.state = init.state, remaining.delays = remaining.delays))
}



Simulation <- function(model, init.state, off.switch, number.simulations) {
'This is one simulation. It takes the model and initial states and for the
chosen number of iterations/timesteps, it performs state transitions and returns the result.
It calls the Intervention functions: Delay, OffSwitch and ResurrectionStop.'
    
  time.steps = 40
  all.transitions <- list()
    
  for (n in 1:number.simulations){
    transitions <- matrix(NA, nrow = time.steps+1, ncol = length(model$nodes))
    colnames(transitions) <- model$nodes
    transitions[1,] <- init.state
    
    # copy delays vector to keep track for remaining delays for each gene
    remaining.delays <- delays   
    off <- off.switch
    timesteps.counter = 0 
    
     for (t in 1:time.steps) {
      trans <- transitions[t, ]
      timesteps.counter <- timesteps.counter + 1
      
      # interventions
      intervention <- Interventions(trans, timesteps.counter, off, remaining.delays, model, init.state)
      remaining.delays <- intervention$remaining.delays
      off <- intervention$off
      
      # perform state transition
      transitions[t + 1, ] <- stateTransition(intervention$updated.model, transitions[t, ], "asynchronous") 
     }
    # list that collects all matrices, access via all.transitions[[time.steps]]
    all.transitions[[n]] <- transitions 
    
    # option to print and plot all transition steps # get rid of that
    #print(all.transitions[[n]])
    #plotSequence(sequence = transitions, onColor = "#4EB7BA", offColor = "#0C2577", borderColor = "black")
    
    # option to plot the path to attractor state
    path = getPathToAttractor(model, init.state)
    plotSequence(sequence = path, onColor = "#4EB7BA", offColor = "#0C2577", borderColor = "black")
  }
  node.names <- colnames(transitions)
  return(list(all.transitions = all.transitions, node.names = node.names))
}


Attractors <- function(all.transitions) {
'This function takes the matrix with all transitions from the Simulation 
function and saves the last row.'
  
  last.rows <- matrix(NA, length(all.transitions), ncol = ncol(all.transitions[[1]]))
  colnames(last.rows) <- colnames(all.transitions[[1]])
  
  for (k in 1:length(all.transitions)) { 
    transitions <- all.transitions[[k]]
    last.rows[k, ] <- transitions[nrow(transitions), ]
  }
  last.rows.data <- data.frame(last.rows)

  # calculate the simulation results that was achieved most often 
  # calculate the percentage of the most frequent row
  num.identical <- max(table(duplicated(last.rows.data)))
  score <- (num.identical+1)*100 / nrow(last.rows.data)
  most.frequent.row <- last.rows.data[which.max(table(last.rows.data)), ]
  
  # option to print the attractor score and the most frequent outcome.
  print(paste("Attractor score: ", score))
  print("Most frequent outcome:")
  print(most.frequent.row)
  
  return(list(last.rows = last.rows.data, score = score, most.frequent.row = most.frequent.row))
}


SenescenceOutcome<- function(most.frequent.row){
'This function checks the senescence node of the attractor and prints out the result.'
  
  if (most.frequent.row[, "sc"] == 0) {
    print("Senescencent cell was eliminated.")
  }
  else {
    print("Senescencent cell was not eliminated.")
  }
}
  

RunModel <- function(model, state) {
'This function takes the model, simulates the process as many times as chosen 
and overlaps the outcomes of each simulation in a plot for visualization.'
  
  number.simulations <- 200
  results.df <- data.frame()
  
  # start simulation and collect results 
  simulation <- Simulation(model, init.state, off.switch, number.simulations)
  all.transitions <- simulation$all.transitions
  node.names <- simulation$node.names
  
  # calculate attractors
  attractor.results <- Attractors(all.transitions)
  
  # evaluate outcome
  score <- attractor.results$score
  most.frequent.row <- attractor.results$most.frequent.row
  evaluate <- SenescenceOutcome(most.frequent.row)
 
 
  # calculate average of each simulation
  average.list <- lapply(seq_len(nrow(all.transitions[[1]])), function(row) {
    lapply(seq_len(ncol(all.transitions[[1]])), function(col) {
      average <- mean(sapply(all.transitions, function(mat) mat[row, col]))
      return(average)
    })
  })
  
  # plot outcome of all simulations
  average.matrix <- matrix(unlist(average.list), nrow = nrow(all.transitions[[1]]), ncol = ncol(all.transitions[[1]]), byrow = TRUE)
  colnames(average.matrix) <-  node.names
  transposed.average <- t(average.matrix)
  valid.values <- transposed.average[!is.na(transposed.average) & is.numeric(transposed.average) & transposed.average >= 0 & transposed.average <= 1]

  heatmap.2(transposed.average, col = colorRampPalette(c("white", "grey", "black"))(100), 
            breaks = seq(min(valid.values), max(valid.values), length.out = 101),
            trace = "none", dendrogram = "none", key = FALSE, Rowv = FALSE, Colv = FALSE,
            margins = c(5, 5), cexRow = 0.8, cexCol = 0.8,
            labRow = node.names)
}

## Run model ##

RunModel(model, state)

