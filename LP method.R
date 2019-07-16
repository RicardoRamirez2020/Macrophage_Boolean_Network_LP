kfoldCV<-function(kfold, times, obs, delta, lambda, Z, b, n, K, T_, annot, annot_node, active_mu, active_sd, 
															 inactive_mu, inactive_sd, mu_type, delta_type, prior=NULL, sourceNode=NULL, sinkNode=NULL,
															 allint=FALSE, allpos=FALSE) {
	
  obs_kfold <- list()
  num <- dim(obs)[1] * dim(obs)[2] * (dim(obs)[3]-1) # observations from t=1 are not removed
  le <- ceiling(num / kfold)
	sq_err_all <- edges_all <- baseline_all <- vector()
	
	# trick so that observations from t=1 are not removed, part 1
	obs_complete <- obs
	obst1 <- obs[ , ,1]
	obs <- obs[ , ,-1]
	
	NAentries <- which(is.na(obs))

	# define k-fold groups: stratified
  for (k in 1:times) {
		tmps <- sample(seq(1, kfold), kfold)
		sq_err_tmp <- vector()
		
		for (j in 1:kfold) {
			tmp <- c(tmps[j:kfold], tmps[1:j-1])
			obs_order <- order(obs)
			obs_kfold[[j]] <- array(NA, c(n,K,T_-1))
			
			for (i in 1:le) {
				if (num >= i)
					obs_kfold[[j]][obs_order[tmp[1]]] <- obs[obs_order[tmp[1]]]

				obs_order <- obs_order[-tmp]
				tmp <- c(tmp[-1], tmp[1])
			}
			obs_kfold[[j]] <- array(obs_kfold[[j]], c(n,K,T_-1))
		}
		
		for (x in 1:kfold) {
			test_ids <- seq(1, kfold)[-x]
			train_tmp <- vector()
			
			for(i in test_ids)
				train_tmp <- rbind(train_tmp, c(obs_kfold[[i]]))
			
			train_data <- rep(NA, dim(train_tmp)[2])
			
			for (i in 1:dim(train_tmp)[2]) {
				if (!all(is.na(train_tmp[,i]))) {
					train_data[i] <- sum(train_tmp[,i], na.rm=T)
				}
			}
			
			train_data <- array(train_data, c(n,K,T_-1))
			obs_modified <- train_data
		
			rem_entries <- which(is.na(obs_modified), arr.ind=TRUE)
			rem_entries_vec <- which(is.na(obs_modified))
			rowsToRemove <- which(rem_entries_vec %in% intersect(rem_entries_vec, NAentries))

			if (length(rowsToRemove) > 0){
				rem_entries <- rem_entries[-rowsToRemove,]
				rem_entries_vec <- rem_entries_vec[-rowsToRemove]
			}
			rem_entries[,3] <- rem_entries[,3] + 1
			rem_entries_vec <- rem_entries_vec + n*K

			# trick so that observations from t=1 are not removed, part 2
			obstemp <- array(NA, c(n,K,T_))
			obstemp[ , ,1] <- obst1
			obstemp[ , ,2:T_] <- obs_modified
			obs_modified <- obstemp

			## do ILP
			res <- ILP(obs=obs_modified, delta=delta, lambda=lambda, Z=Z, b=b, n=n, K=K, T_=T_, annot, delta_type=delta_type,
															 prior=prior, sourceNode=sourceNode,sinkNode=sinkNode,all.int=allint,all.pos=allpos)

			adja <- getAdja(res=res, n=n)
			baseline <- getBaseline(res=res, n=n)
			edges_all <- rbind(edges_all, c(t(adja)))
			baseline_all <- rbind(baseline_all, baseline)
			
			# predict missing observation values
			predict <- .calcPredictionKfoldCV_timeSeries(obs=obs_modified, delta=delta, b=b, n=n, K=K, adja=adja, baseline=baseline, 
																									 rem_entries=rem_entries, rem_entries_vec=rem_entries_vec, active_mu=active_mu,
																									 active_sd=active_sd, inactive_mu=inactive_mu, inactive_sd=inactive_sd, mu_type=mu_type)

			ids_rem <- rem_entries_vec
			sq_err_tmp <- c(sq_err_tmp, ((predict[ids_rem] - obs_complete[ids_rem])^2))  # calculate squared error of predicted and observed
		}
		sq_err_all <- rbind(sq_err_all, sq_err_tmp)
  } # end times
  
  sq_err <- apply(sq_err_all, 2, mean, na.rm=T)
  tmp1 <- rep(annot_node, rep(n,n))
  tmp2 <- rep(annot_node, n)
  id_selfloop <- which(tmp1 == tmp2)
  tmp <- paste(tmp1, tmp2, sep="->")
  edges_all <- edges_all[ ,-id_selfloop]
  colnames(edges_all) <- tmp[-id_selfloop]
  MSE <- mean(sq_err, na.rm=T)
  
  return(list(MSE=MSE, edges_all=edges_all, baseline_all=baseline_all))
}

.calcPredictionKfoldCV_timeSeries <- function(obs, delta, b ,n ,K, adja, baseline, rem_entries, rem_entries_vec,
																						  active_mu, active_sd, inactive_mu, inactive_sd, mu_type) {

	inact_entries <- which(b == 0)
	predict <- obs

	if (dim(rem_entries)[1] > 0) {
		if (mu_type == "single") {
		
			for (ent in 1:dim(rem_entries)[1]) {
				rem_gene <- rem_entries[ent,1]
				rem_k <- rem_entries[ent,2]
				rem_t <- rem_entries[ent,3]
				
				rem_ent_test <- rem_entries_vec[ent] %% (n*K)
				if (rem_ent_test == 0) 
					rem_ent_test <- n*K  # when the removed entry is in the last position the remainder of the division is zero and the test fails
				res <- rem_ent_test %in% inact_entries
				
				delta_rem <- delta[rem_gene]
				predict <- .calculatePredictionValue_Kfold_ts(predict, obs, n, adja, baseline, delta, delta_rem, b, active_mu, active_sd, 
																											inactive_mu, inactive_sd, rem_gene, rem_k, rem_t, res, mu_type)
			}
		}
		else if (mu_type == "perGene") {

			for (ent in 1:dim(rem_entries)[1]) {
				rem_gene <- rem_entries[ent,1]
				rem_k <- rem_entries[ent,2]
				rem_t <- rem_entries[ent,3]
				
				rem_ent_test <- rem_entries_vec[ent] %% (n*K)
				if (rem_ent_test == 0) 
					rem_ent_test <- n*K  # when the removed entry is in the last position the remainder of the division is zero and the test fails
				res <- rem_ent_test %in% inact_entries
				
				delta_rem <- delta[rem_gene]
				predict <- .calculatePredictionValue_Kfold_ts(predict, obs, n, adja, baseline, delta, delta_rem, b, active_mu[rem_gene], 
																											active_sd[rem_gene], inactive_mu[rem_gene], inactive_sd[rem_gene], 
																											rem_gene, rem_k, rem_t, res, mu_type)
			}
		}
		else if (mu_type == "perGeneExp") {
		
			for (ent in 1:dim(rem_entries)[1]) {
				rem_gene <- rem_entries[ent,1]
				rem_k <- rem_entries[ent,2]
				rem_t <-rem_entries[ent,3]
				
				rem_ent_test = rem_entries_vec[ent] %% (n*K)
				if (rem_ent_test == 0) 
					rem_ent_test <- n*K # when the removed entry is in the last position the remainder of the division is zero and the test fails
				res <- rem_ent_test %in% inact_entries
				
				delta_rem <- delta[rem_gene, rem_k]
				predict <- .calculatePredictionValue_Kfold_ts(predict, obs, n, adja, baseline, delta, delta_rem, b, active_mu[rem_gene, rem_k], 
																											active_sd[rem_gene, rem_k], inactive_mu[rem_gene, rem_k], inactive_sd[rem_gene, rem_k], 
																											rem_gene, rem_k, rem_t, res, mu_type)
			}
		}
		
		else if (mu_type == "perGeneTime") {
		
			for (ent in 1:dim(rem_entries)[1]) {
				rem_gene <- rem_entries[ent,1]
				rem_k <- rem_entries[ent,2]
				rem_t <- rem_entries[ent,3]
				
				rem_ent_test <- rem_entries_vec[ent]%%(n*K)
				if (rem_ent_test == 0) 
					rem_ent_test <- n*K  # when the removed entry is in the last position the remainder of the division is zero and the test fails
				res <- (rem_ent_test %in% inact_entries)
				
				delta_rem <- delta[rem_gene, rem_t]
				predict <- .calculatePredictionValue_Kfold_ts(predict, obs, n, adja, baseline, delta, delta_rem, b,
																											active_mu[rem_gene, rem_t], active_sd[rem_gene, rem_t], 
																											inactive_mu[rem_gene, rem_t], inactive_sd[rem_gene, rem_t], 
																											rem_gene, rem_k, rem_t, res, mu_type)
			}
		}
		
		else if (mu_type == "perGeneExpTime") {
		
			for (ent in 1:dim(rem_entries)[1]) {
				rem_gene <- rem_entries[ent,1]
				rem_k <- rem_entries[ent,2]
				rem_t <- rem_entries[ent,3]
				
				rem_ent_test = rem_entries_vec[ent]%%(n*K)
				if (rem_ent_test == 0) 
					rem_ent_test <- n*K  # when the removed entry is in the last position the remainder of the division is zero and the test fails
				res <- rem_ent_test %in% inact_entries
				
				delta_rem <- delta[rem_gene, rem_k, rem_t]
				predict <- .calculatePredictionValue_Kfold_ts(predict, obs, n, adja, baseline, delta, delta_rem, b,
																											active_mu[rem_gene, rem_k, rem_t], active_sd[rem_gene, rem_k, rem_t], 
																											inactive_mu[rem_gene, rem_k, rem_t], inactive_sd[rem_gene, rem_k, rem_t], 
																											rem_gene, rem_k, rem_t, res, mu_type)
			}
		}
		
	}
	return(predict)
}


.calculatePredictionValue_Kfold_ts <- function(predict, obs, n, adja, baseline, delta, delta_rem, b, active_mu, active_sd, 
																							 inactive_mu, inactive_sd, rem_gene, rem_k, rem_t, res, mu_type) {

	# if the removed entry is an inactive node due to some knockdown, then predict as inactivet=rem_entries[ent,3]
	if (res == TRUE) {
		predict[rem_gene,rem_k,rem_t] <- rnorm(1, inactive_mu, inactive_sd)
	}
	else {
		if (is.na(baseline[rem_gene])) {
			in_flow <- 0
		}
		else {
			in_flow <- baseline[rem_gene]
		}
		
		pa <- which(adja[ ,rem_gene] != 0)
		if (length(pa) == 0) {  # if there are no parents: rem_gene is root node 
			if (in_flow >= delta_rem) {  # root node is active if its inflow is greater than its delta
				predict[rem_gene,rem_k,rem_t] <- rnorm(1, active_mu, active_sd) 
			}
			else {
				predict[rem_gene,rem_k,rem_t] <- rnorm(1, inactive_mu, inactive_sd)
			}
		}
		else {
			flagNA <- 0
			for (j in 1:length(pa)) {
				
				delta_pa <- .setDeltaValue_KfoldCV_ts(delta, pa, j, rem_k, rem_t, mu_type)
				if (is.na(obs[pa[j],rem_k,rem_t-1])) {  # if parent observation is NA
					flagNA <- 1
					
					if ((adja[pa[j],rem_gene] < 0) & (b[(rem_k-1)*n + pa[j]] == 1)) {  # if the incoming edge is negative and the parent node is active, node state is unknown
						predict[rem_gene,rem_k,rem_t] <- NA
						return(predict=predict)
					}
				}
				else if ((obs[pa[j],rem_k,rem_t-1] >= delta_pa) & (b[(rem_k-1)*n + pa[j]] == 1)) {  # if parent is active and not silenced calculate node inflow
					in_flow <- sum(in_flow, adja[pa[j],rem_gene] * obs[pa[j],rem_k,rem_t-1], na.rm=T)
				}
			}
			
			if ((flagNA == 1) & (in_flow < delta_rem)) {  # if parent observation is NA and inflow < delta, node state is unknown
				predict[rem_gene,rem_k,rem_t] <- NA
				return(predict)
			}			
			
			# predict node state according to inflow
			if (in_flow >= delta_rem) {  
				predict[rem_gene,rem_k,rem_t] <- rnorm(1, active_mu, active_sd) 
			}
			else {
				predict[rem_gene,rem_k,rem_t] <- rnorm(1, inactive_mu, inactive_sd)
			}
		}
	}
	return(predict)
}

.setDeltaValue_KfoldCV_ts <- function(delta, pa, j, rem_k, rem_t, mu_type) {

	if ((mu_type == "single") | (mu_type == "perGene")) {
		delta_pa <- delta[pa[j]]
	}
	else if (mu_type == "perGeneExp") {
		delta_pa <- delta[pa[j], rem_k]
	}
	else if (mu_type == "perGeneTime") {
		delta_pa <- delta[pa[j], rem_t-1]
	}
	else if (mu_type == "perGeneExpTime") {
		delta_pa <- delta[pa[j], rem_k, rem_t-1]
	}

	return(delta_pa)
}

ILP <- function(obs, delta, lambda, Z, b, n, K, T_, annot, delta_type, prior=NULL, sourceNode=NULL,
															sinkNode=NULL, all.int=FALSE, all.pos=FALSE) {

	nConstr <- n*K*(T_-1)
	
  ## weight matrix of dim ((K*n)x(2n²+n)) (w_i^0)
  if(all.pos) {
		W <- matrix(0, nrow=nConstr, ncol=n*n+n)
	}
  else {
		W <- matrix(0, nrow=nConstr, ncol=2*n*n+n)
	}
  colnames(W) <- annot
  
  # direction of inequation
  f.dir <- rep("<=", nConstr)
  
  # Vector of numeric values for the right-hand sides of the constraints
  bvec <- rep(0, nConstr)
  J <- seq(1,n)
  count <- 1
  
  if(delta_type == "perGene"){
  
		if (all.int)
			delta <- rep(1, n)
		
		for (t in 2:T_) {
			for (k in 1:K) {
				for (i in 1:n) {
					
					if (b[(k-1)*n + i] == 1) {  # if the entry in b is 1 then the gene is active
						if (!is.na(obs[i,k,t])) {  # if the observation=NA, just do nothing
							
							if (obs[i,k,t] >= delta[i]) {  # if observation of gene i after knockdown k is active
								if (all.pos) {
									W[count,i+(n*n)] <- 1  # set offset parameter (baseline of gene i)
									
									for (j in J[J!=i]) {  # sum
										idPos <- which(annot == paste("w+", j, i, sep="_"))  # positive parameter
										
										if (!is.na(obs[j,k,t-1])) {
											if (obs[j,k,t-1] >= delta[j] & b[(k-1)*n + j] == 1) {
												W[count,idPos] <- obs[j,k,t-1]
											}
											else {
												W[count,idPos] <- 0
											}
										}
										else {
											W[count,idPos] <- NA
										}
									}
								}
								
								else {
									W[count,i+(2*n*n)] <- 1
									
									for (j in J[J!=i]) {  # sum
										idPos <- which(annot == paste("w+", j, i, sep="_"))  # positive parameter
										idNeg <- which(annot == paste("w-", j, i, sep="_"))  # negative parameter
										
										if (!is.na(obs[j,k,t-1])) {
											if (obs[j,k,t-1] >= delta[j] & b[(k-1)*n + j] == 1) {
												W[count,idPos] <- obs[j,k,t-1]
												W[count,idNeg] <- -obs[j,k,t-1]
											}
											else {
												W[count,idPos] <- 0
												W[count,idNeg] <- 0
											}
										}
										else {
											W[count,idPos] <- NA
											W[count,idNeg] <- NA
										}
									} 
								}
								f.dir[count] <- ">="
								bvec[count] <- delta[i]
							}
							
							if (obs[i,k,t]< delta[i]) {  # if observation of gene i after knockdown k is NOT acitve
								if (all.pos) {
									W[count,i+(n*n)] <- 1  # set offset parameter (baseline of gene i)
									
									for (j in J[J!=i]) {  # sum
										idPos <- which(annot == paste("w+", j, i, sep="_"))  # positive parameter
										
										if (!is.na(obs[j,k,t-1])) {
											if (obs[j,k,t-1] >= delta[j] & b[(k-1)*n + j] == 1) {
												W[count,idPos] <- obs[j,k,t-1]
											}
											else {
												W[count,idPos] <- 0
											}
										}
										else {
											W[count,idPos] <- NA
										}
									}
								}
								else {
									W[count,i+(2*n*n)] <- 1
									
									for (j in J[J!=i]) {  # sum
										idPos <- which(annot == paste("w+", j, i, sep="_"))  # positive parameter
										idNeg <- which(annot == paste("w-", j, i, sep="_"))  # negative parameter
										
										if (!is.na(obs[j,k,t-1])) {
											if (obs[j,k,t-1] >= delta[j] & b[(k-1)*n + j] == 1) {
												W[count,idPos] <- obs[j,k,t-1]
												W[count,idNeg] <- -obs[j,k,t-1]
											}
											else {
												W[count,idPos] <- 0
												W[count,idNeg] <- 0
											}
										}
										else {
											W[count,idPos] <- NA
											W[count,idNeg] <- NA
										}
									}
								}
							f.dir[count] <- "<="
							bvec[count] <- 0
							}
						}
					}
					count <- count+1
				} # end i
			} # end k
		} # end t
	}
	else if (delta_type == "perGeneExp") {
  
		if(all.int)
			delta <- matrix(rep(1,n*K), nrow=n, ncol=K)
			
		for (t in 2:T_) {
			for (k in 1:K) {
				for (i in 1:n) {
					
					if (b[(k-1)*n + i] == 1) {  # if the entry in b is 1 then the gene is active
						if (!is.na(obs[i,k,t])) {  # if the observation=NA, just do nothing
							
							if (obs[i,k,t] >= delta[i,k]) {  # if observation of gene i after knockdown k is active
								if (all.pos) {
									W[count,i+(n*n)] <- 1  # set offset parameter (baseline of gene i)
									
									for (j in J[J!=i]) {  # sum
										idPos <- which(annot == paste("w+", j, i, sep="_"))  # positive parameter
										
										if (!is.na(obs[j,k,t-1])) {
											if (obs[j,k,t-1] >= delta[j,k] & b[(k-1)*n + j] == 1) {
												W[count,idPos] <- obs[j,k,t-1]
											}
											else {
												W[count,idPos] <- 0
											}
										}
										else {
											W[count,idPos] <- NA
										}
									}
								}
								else {
									W[count,i+(2*n*n)] <- 1
									for (j in J[J!=i]) {
										idPos <- which(annot == paste("w+", j, i, sep="_"))
										idNeg <- which(annot == paste("w-", j, i, sep="_"))
										
										if (!is.na(obs[j,k,t-1])) {
											if (obs[j,k,t-1] >= delta[j,k] & b[(k-1)*n + j] == 1) {
												W[count,idPos] <- obs[j,k,t-1]
												W[count,idNeg] <- -obs[j,k,t-1]
											}
											else {
												W[count,idPos] <- 0
												W[count,idNeg] <- 0
											}
										}
										else {
											W[count,idPos] <- NA
											W[count,idNeg] <- NA
										}
									} 
								}
								f.dir[count] <- ">="
								bvec[count] <- delta[i,k]
							}
							
							if (obs[i,k,t]< delta[i,k]) {  # if observation of gene i after knockdown k is NOT acitve
								if (all.pos) {  # set offset parameter (baseline of gene i)
									W[count,i+(n*n)] <- 1
									
									for (j in J[J!=i]) {  # sum
										idPos <- which(annot == paste("w+", j, i, sep="_"))  # positive parameter
										
										if (!is.na(obs[j,k,t-1])) {
											if (obs[j,k,t-1] >= delta[j,k] & b[(k-1)*n + j] == 1){
												W[count,idPos] <- obs[j,k,t-1]
											}
											else {
												W[count,idPos] <- 0
											}
										}
										else {
											W[count,idPos] <- NA
										}
									}
								}
								else {
									W[count,i+(2*n*n)] <- 1
									
									for (j in J[J!=i]) {  # sum
										idPos <- which(annot == paste("w+", j, i, sep="_"))  # positive parameter
										idNeg <- which(annot == paste("w-", j, i, sep="_"))  # negative parameter
										
										if (!is.na(obs[j,k,t-1])) {
											if (obs[j,k,t-1] >= delta[j,k] & b[(k-1)*n + j] == 1) {
												W[count,idPos] <- obs[j,k,t-1]
												W[count,idNeg] <- -obs[j,k,t-1]
											}
											else {
												W[count,idPos] <- 0
												W[count,idNeg] <- 0
											}
										}
										else {
											W[count,idPos] <- NA
											W[count,idNeg] <- NA
										}
									}
								}
							f.dir[count] <- "<="
							bvec[count] <- 0
							}
						}
					}
					count <- count + 1
				} # end i
			} # end k
		} # end t
	}
	else if (delta_type == "perGeneTime"){
  
		if (all.int)
			delta <- matrix(rep(1,n*T_), nrow=n, ncol=T_)
			
		for (t in 2:T_) {
			for (k in 1:K) {
				for (i in 1:n) {
					
					if (b[(k-1)*n + i] == 1) {  # if the entry in b is 1 then the gene is active
						if (!is.na(obs[i,k,t])) {  # if the observation=NA, just do nothing
							
							if (obs[i,k,t] >= delta[i,t]) {  # if observation of gene i after knockdown k is active
								if (all.pos) {
									W[count,i+(n*n)] <- 1  # set offset parameter (baseline of gene i)
									
									for (j in J[J!=i]) {  # sum
										idPos <- which(annot == paste("w+", j, i, sep="_"))  # positive parameter
										
										if (!is.na(obs[j,k,t-1])) {
											if (obs[j,k,t-1] >= delta[j,t-1] & b[(k-1)*n + j] == 1) {
												W[count,idPos] <- obs[j,k,t-1]
											}
											else {
												W[count,idPos] <- 0
											}
										}
										else {
											W[count,idPos] <- NA
										}
									}
								}
								else {
									W[count,i+(2*n*n)] <- 1
									for (j in J[J!=i]) {
										idPos <- which(annot == paste("w+", j, i, sep="_"))
										idNeg <- which(annot == paste("w-", j, i, sep="_"))
										
										if (!is.na(obs[j,k,t-1])) {
											if(obs[j,k,t-1] >= delta[j,t-1] & b[(k-1)*n + j] == 1){
												W[count,idPos] <- obs[j,k,t-1]
												W[count,idNeg] <- -obs[j,k,t-1]
											}
											else {
												W[count,idPos] <- 0
												W[count,idNeg] <- 0
											}
										}
										else {
											W[count,idPos] <- NA
											W[count,idNeg] <- NA
										}
									} 
								}
								f.dir[count] <- ">="
								bvec[count] <- delta[i,t]
							}
							
							if (obs[i,k,t] < delta[i,t]) {  # if observation of gene i after knockdown k is NOT acitve
								if (all.pos) {
									W[count,i+(n*n)] <- 1  # set offset parameter (baseline of gene i)
									
									for (j in J[J!=i]) {  # sum
										idPos <- which(annot == paste("w+", j, i, sep="_"))  # positive parameter
										
										if (!is.na(obs[j,k,t-1])) {
											if (obs[j,k,t-1] >= delta[j,t-1] & b[(k-1)*n + j] == 1){
												W[count,idPos] <- obs[j,k,t-1]
											}
											else {
												W[count,idPos] <- 0
											}
										}
										else {
											W[count,idPos] <- NA
										}
									}
								}
								else {
									W[count,i+(2*n*n)] <- 1
									
									for (j in J[J!=i]) {  # sum
										idPos <- which(annot == paste("w+", j, i, sep="_"))  # positive parameter
										idNeg <- which(annot == paste("w-", j, i, sep="_"))  # negative parameter
										
										if (!is.na(obs[j,k,t-1])) {
											if (obs[j,k,t-1] >= delta[j,t-1] & b[(k-1)*n + j] == 1) {
												W[count,idPos] <- obs[j,k,t-1]
												W[count,idNeg] <- -obs[j,k,t-1]
											}
											else {
												W[count,idPos] <- 0
												W[count,idNeg] <- 0
											}
										}
										else {
											W[count,idPos] <- NA
											W[count,idNeg] <- NA
										}
									}
								}
							f.dir[count] <- "<="
							bvec[count] <- 0
							}
						}
					}
					count <- count+1
				} # end i
			} # end k
		} # end t
	}
	else if (delta_type == "perGeneExpTime") {
  
		if (all.int)
			delta <- array(rep(1,n*K*T_), c(n,K,T_))
		
		for (t in 2:T_) {
			for (k in 1:K) {
				for (i in 1:n) {
					
					if (b[(k-1)*n + i] == 1) {  # if the entry in b is 1 then the gene is active
						if (!is.na(obs[i,k,t])) {  # if the observation=NA, just do nothing
							
							if (obs[i,k,t] >= delta[i,k,t]) {  # if observation of gene i after knockdown k is active
								if (all.pos) {
									W[count,i+(n*n)] <- 1  # set offset parameter (baseline of gene i)
									
									for (j in J[J!=i]) {  # sum
										idPos <- which(annot == paste("w+", j, i, sep="_"))  # positive parameter
										
										if (!is.na(obs[j,k,t-1])) {
											if (obs[j,k,t-1] >= delta[j,k,t-1] & b[(k-1)*n + j] == 1) {
												W[count,idPos] <- obs[j,k,t-1]
											}
											else {
												W[count,idPos] <- 0
											}
										}
										else {
											W[count,idPos] <- NA
										}
									}
								}
								else {
									W[count,i+(2*n*n)] <- 1
									for (j in J[J!=i]) {
										idPos <- which(annot == paste("w+", j, i, sep="_"))  # positive parameter
										idNeg <- which(annot == paste("w-", j, i, sep="_"))  # negative parameter
										
										if (!is.na(obs[j,k,t-1])) {
											if (obs[j,k,t-1] >= delta[j,k,t-1] & b[(k-1)*n + j] == 1) {
												W[count,idPos] <- obs[j,k,t-1]
												W[count,idNeg] <- -obs[j,k,t-1]
											}
											else {
												W[count,idPos] <- 0
												W[count,idNeg] <- 0
											}
										}
										else {
											W[count,idPos] <- NA
											W[count,idNeg] <- NA
										}
									} 
								}
								f.dir[count] <- ">="
								bvec[count] <- delta[i,k,t]
							}
							
							if (obs[i,k,t] < delta[i,k,t]) {  # if observation of gene i after knockdown k is NOT acitve
								if (all.pos) {
									W[count,i+(n*n)] <- 1  # set offset parameter (baseline of gene i)
									
									for (j in J[J!=i]) {  # sum
										idPos <- which(annot == paste("w+", j, i, sep="_"))  # positive parameter
										
										if (!is.na(obs[j,k,t-1])) {
											if (obs[j,k,t-1] >= delta[j,k,t-1] & b[(k-1)*n + j] == 1) {
												W[count,idPos] <- obs[j,k,t-1]
											}
											else {
												W[count,idPos] <- 0
											}
										}
										else {
											W[count,idPos] <- NA
										}
									}
								}
								else {
									W[count,i+(2*n*n)] <- 1
									
									for (j in J[J!=i]) {  # sum
										idPos <- which(annot==paste("w+",j,i,sep="_"))  # positive parameter
										idNeg <- which(annot==paste("w-",j,i,sep="_"))  # negative parameter
										
										if (!is.na(obs[j,k,t-1])) {
											if (obs[j,k,t-1] >= delta[j,k,t-1] & b[(k-1)*n + j] == 1) {
												W[count,idPos] <- obs[j,k,t-1]
												W[count,idNeg] <- -obs[j,k,t-1]
											}
											else {
												W[count,idPos] <- 0
												W[count,idNeg] <- 0
											}
										}
										else {
											W[count,idPos] <- NA
											W[count,idNeg] <- NA
										}
									}
								}
							f.dir[count] <- "<="
							bvec[count] <- 0
							}
						}
					}
					count <- count+1
				} # end i
			} # end k
		} # end t
	}

  
  # now add slack varibles to W: 
  if (lambda != 0) {
		sl <- matrix(0, nrow=nConstr, ncol=nConstr)
		annot_s <- paste("s",seq(1, nConstr), sep="_")
		colnames(sl) <- annot_s
		
		# attention: for the constr. where observation is smaller than threshold
		xi <- vector()
		for (j in 1:length(f.dir)) {
			if(f.dir[j] == ">=") 
				xi <- c(xi ,0)
			if(f.dir[j] == "<=") 
				xi <- c(xi, -1)
		}
		diag(sl) <- xi
		W <- cbind(W, sl)
		if (all.pos) {
			cvec <- c(rep(1, n*n), rep(1/Z, n), rep(1/lambda, nConstr))
			# self-activation is not allowed
			id_self <- c(which(annot == paste("w+", seq(1, n), seq(1, n), sep="_")))
		}
		else {
			cvec <- c(rep(1, n*n), rep(1, n*n), rep(1/Z, n), rep(1/lambda, nConstr))
			# self-activation is not allowed
			id_self <- c(which(annot == paste("w+", seq(1, n) ,seq(1, n), sep="_")),
									 which(annot == paste("w-", seq(1, n), seq(1, n), sep="_")))
		}
		cvec[id_self] <- 0
		names(cvec) <- c(annot, annot_s)
		
	}
	else {
		if (all.pos) {
			cvec <- c(rep(1, n*n), rep(1/Z, n)) 
			# self-activation is not allowed
			id_self <- c(which(annot == paste("w+", seq(1, n), seq(1, n), sep="_")))
		}
		else {
			cvec <- c(rep(1, n*n), rep(1, n*n), rep(1/Z, n)) 
			# self-activation is not allowed
			id_self <- c(which(annot == paste("w+", seq(1, n), seq(1, n), sep="_")),
									 which(annot == paste("w-", seq(1, n), seq(1, n), sep="_")))
		}
		cvec[id_self] <- 0
		names(cvec) <- c(annot)
  }

  # condition that each node which is not End hast at least delta[i] outgoing edges
  if (!is.null(sinkNode)) {
		W_tmp1 <- vector()
		gene_tmp <- seq(1, n)[-sinkNode]
		
		for (i in gene_tmp) {
			# outgoing edge can come from all nodes except itself
			tmp <- seq(1, n)[-i]
			
			if (length(tmp) > 1) {
				# for negative and positive parameter
				annot_pos <- paste("w+", i, tmp, sep="_")
				
				if (!all.pos) 
					annot_neg <- paste("w-", i, tmp, sep="_")
				
				add_row <- rep(0, length(cvec))
				add_row[which(annot %in% annot_pos)] <- 1
				
				if (!all.pos) 
					add_row[which(annot %in% annot_neg)] <- 1
				
				W_tmp1 <- rbind(W_tmp1, as.double(add_row))
				bvec <- c(bvec, delta[i])
				f.dir <- c(f.dir, ">=")
			}
		}
		W <- rbind(W, W_tmp1)
  }
  
  # conditions that each node which is not Start has at least delta[i] incoming edges
  if (!is.null(sourceNode)) {
		W_tmp2 <- vector()
		gene_tmp <- seq(1, n)[-sourceNode]
		
		for (i in gene_tmp) {
			# incoming edge can come from all nodes except itself
			tmp <- seq(1, n)[-i]
			
			if (length(tmp) > 1) {
				annot_pos <- paste("w+", tmp, i, sep="_")
				if (!all.pos) 
					annot_neg <- paste("w-", tmp, i, sep="_")
				
				add_row <- rep(0, length(cvec))
				add_row[which(annot %in% annot_pos)] <- 1
				
				if (!all.pos) 
					add_row[which(annot %in% annot_neg)] <- 1
				
				W_tmp2 <- rbind(W_tmp2, as.double(add_row))
				bvec <- c(bvec, delta[i])
				f.dir <- c(f.dir, ">=")
			}
		}
		W <- rbind(W, W_tmp2)
  }

  ## if there is a prior
  if (!is.null(prior)) {
		for (i in 1:length(prior)) {
		
			tmp <- rep(0, dim(W)[2])
			tmp[which(prior[[i]][1] == annot)] <- as.double(prior[[i]][2])
			W <- rbind(W, tmp)
			bvec <- c(bvec, as.double(prior[[i]][4]))
			f.dir <- c(f.dir, prior[[i]][3])
		}
  }

  # Maximize the gross margin
  # min - direction of optimization
  # cvec - objective function (Numeric vector of coefficients of objective function)
  # W - Matrix of numeric constraint coefficients, one row per constraint, one column per variable
  # f.dir vector of character strings giving the direction of the constraint
  # bvec - vector of numeric values for the right-hand sides of the constraints
  res <- lp("min", cvec, W, f.dir, bvec, all.int=all.int) 

  return(res)
}


library('lpNet')
library("R.matlab")

n <- 27
K <- 4
T_ <- 50
delta_type <- "perGeneExp"

LPS <- readMat("C:/Users/RJ/Documents/college/manuscript/LPS/LPS.mat")
LPS <- array(unlist(LPS), dim = c(27, 4, 50))
B <- readMat("C:/Users/RJ/Documents/college/manuscript/B.mat")
B <- array( unlist(B), dim = c(27, 4, 50))
LPSmean <- readMat("C:/Users/RJ/Documents/college/manuscript/LPS/LPSmean.mat")
LPSmean <- array( unlist(LPSmean), dim = c(27, 4))
names <- readMat("C:/Users/RJ/Documents/college/manuscript/names.mat")
names <- c( unlist(names))

active_mu <- readMat("C:/Users/RJ/Documents/college/manuscript/LPS/mu_active.mat")
active_mu <- array(unlist(active_mu), dim = c(27, 4))
active_sd <- readMat("C:/Users/RJ/Documents/college/manuscript/LPS/sd_active.mat")
active_sd <- array(unlist(active_sd), dim = c(27, 4))
inactive_mu <- readMat("C:/Users/RJ/Documents/college/manuscript/LPS/mu_inactive.mat")
inactive_mu <- array(unlist(inactive_mu), dim = c(27, 4))
inactive_sd <- readMat("C:/Users/RJ/Documents/college/manuscript/LPS/sd_inactive.mat")
inactive_sd <- array(unlist(inactive_sd), dim = c(27, 4))

mu_type = "perGeneExp"

times <- 5 
annot_node <- seq(1,n)
lambda <- calcRangeLambda(LPS, LPSmean, delta_type,flag_time_series=TRUE)
Q = lambda
MSE <- Inf
kfold <- 5


for (lamd in lambda) {
kcv_res <- kfoldCV(kfold, times, LPS, LPSmean, lambda=lamd,Z = 100,B, n, K, T_, annot=getEdgeAnnot(n),annot_node, active_mu, active_sd,inactive_mu, inactive_sd, mu_type,delta_type,)
if (kcv_res$MSE < MSE) {
MSE <- kcv_res$MSE
edges_all <- kcv_res$edges_all
bestLambda <- lamd
}
}

res1 <- ILP(LPS, LPSmean, lambda=100, Z = 100,B, n, K, T_ , annot = getEdgeAnnot(n), delta_type)
adja1 <- getAdja(res1, n)
getBaseline(res1,n)
res1
adja1


write.table(adja1, "C:/Users/RJ/Documents/college/manuscript/LPS/adja100.txt", sep="\t")
write.table(getBaseline(res1,n), "C:/Users/RJ/Documents/college/manuscript/LPS/Base100.txt", sep="\t")




