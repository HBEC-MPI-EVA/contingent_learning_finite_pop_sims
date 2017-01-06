# simple cultural accumulation by contingent learning
# this version has flexible group size but density dependent mortality
# this allows flexible migration structure

# my sample function, which doesn't do that dumb integer thing
ersample <- function(x,...) {
    if ( length(x)==1 )
        return(x)
    else {
        return(sample(x,...))
    }
}

# unique that excludes NA
erunique <- function(x) {
    y <- unique(x)
    y[!is.na(y)]
}

sim_accum_flex <- function( tmax=1000 , N=40 , N_init=N/2 , N_groups=2 , 
    s=0.1 , f=0.3, d=1 , e=0 , m=0.01 , death=0.1 , K=0.02 , b=0.1 , u=0 , v=1 , 
    n_behaviors=1e4 , n_teachers=1 ) {
    # tmax: number of time steps
    # N: initial group size
    # N_groups: number of groups of size N
    # n_behaviors: number of unique possible adaptive behaviors
    # s: probability innovation produces a successful solution
    # f: probability of successfully copying a teacher (transmission fidelity)
    # d: probability social learner rejects maladaptive behavior
    # e: probability social learner rejects adaptive behavior
    # m: migration rate
    # u: prob environment changes (independent for each group)
    # v: prob behavior still valuable in new group
    # death: prob any individual dies in each time step
    # r: growth rate
    # n_teachers: number of teachers to sample; if any has non-zero behavior, copy

    # make initial population matrix
    # dimensions are: time t , group g , individual i , feature j (behavior,age)
    # value '0' means any non-adaptive behavior
    pop <- array( NA , dim=c(tmax,N_groups,N,2) )
    # assign behavior to 0
    pop[1,,1:N_init,1] <- 0
    # assign ages to 1
    pop[1,,1:N_init,2] <- 1

    # loop generations
    for ( t in 2:tmax ) {

        # learn
        for ( g in 1:N_groups ) {
            for ( i in 1:N ) {
                if ( !is.na(pop[t-1,g,i,1]) ) {
                    if ( pop[t-1,g,i,1]==0 ) {
                        # do some learning
                        # sample teacher
                        available_teachers <- which(!is.na(pop[t-1,g,,1]))
                        available_teachers <- available_teachers[available_teachers!=i]
                        j <- ersample( available_teachers , size=min(n_teachers,length(available_teachers)) )
                        # check if any have non-zero
                        if ( any(pop[t-1,g,j,1] > 0) ) {
                            # copy one of the non-zero teachers
                            idx <- which(pop[t-1,g,j,1] > 0)
                            j <- ersample( j[idx] , size=1 )
                            pop[t,g,i,1] <- pop[t-1,g,j,1]
                        } else {
                            # no successful teacher
                            pop[t,g,i,1] <- 0
                        }
                        
            			# determine fidelity of transmission
                        #print(f)
            			if ( runif(1) > f ) {
            			     #behavior not successfully copied
            			     pop[t,g,i,1] <- 0
            			}

                        # test if successful
                        do_innovate <- FALSE
                        if ( pop[t,g,i,1] == 0 ) {
                            if ( runif(1) < d ) {
                                # reject unsuccessful behavior
                                do_innovate <- TRUE
                            }
                        } else {
                            if ( runif(1) < e ) {
                                # reject successful behavior
                                do_innovate <- TRUE
                            }
                        }
                        if ( do_innovate==TRUE ) {
                            if ( runif(1) < s ) {
                                # successful innovation
                                pop[t,g,i,1] <- ersample( 1:n_behaviors , size=1 )
                            } else {
                                # unsuccessful innovation
                                pop[t,g,i,1] <- 0
                            }
                        }
                    } else {
                        # behavior not zero, so persist to next time period
                        pop[t,g,i,1] <- pop[t-1,g,i,1]
                    }
                }#!is.na
            }#i
        }#g

        # mortality and aging
        for ( g in 1:N_groups ) {
            # local pop regulation
            # compute mortality prob in group g
            Ng <- sum(!is.na(pop[t,g,,1]))
            pr_d <- min( death + exp(K*Ng)-1 , 1 )
            for ( i in 1:N ) {
                if ( !is.na(pop[t,g,i,1]) ) {
                    if ( runif(1) < pr_d ) {
                        pop[t,g,i,1] <- NA
                        pop[t,g,i,2] <- NA
                    } else {
                        # update age
                        pop[t,g,i,2] <- pop[t-1,g,i,2] + 1
                    }
                    # give birth?
                    if ( any(is.na(pop[t,g,,1])) ) {
                        if ( runif(1) < b ) {
                            j <- ersample( which(is.na(pop[t,g,,1])) , size=1 )
                            pop[t,g,j,1] <- 0
                            pop[t,g,j,2] <- 1
                        }
                    }#birth
                }#!is.na
            }#i
        }#g

        # environmental stochasticity
        # all groups in same ecology, otherwise have to keep track of a lot of stuff
        if ( FALSE ) {
            # old stochasticity code
            for ( g in 1:N_groups ) {
                if ( runif(1) < u ) {
                    # change in env for group g, so test each behavior
                    unique_behaviors <- unique( pop[t,g,,1] )
                    for ( bb in unique_behaviors ) {
                        if ( TRUE ) {
                            # set all instances of b to 0
                            idx <- which( pop[t,g,,1]==bb )
                            if ( length(idx)>0 )
                                pop[t,g,idx,1] <- 0
                        }#v
                    }#b
                }#u
            }#g
        } else {
            # new stochasticity code
            if ( runif(1) < u ) {
                for ( g in 1:N_groups ) {
                    idx <- which( !is.na(pop[t,g,,1]) )
                    if ( length(idx)>0 ) pop[t,g,idx,1] <- 0
                }
            }
        }

        # migrate
        if ( m > 0 & any(is.na(pop[t,,,1])) ) {
            for ( g in 1:N_groups ) {
                for ( i in 1:N ) {
                    if (!is.na(pop[t,g,i,1])) {
                        # check for migration
                        if ( runif(1) < m & any(is.na(pop[t,-g,,1])) ) {
                            ff <- pop[t,g,i,]
                            # find open slot in another group
                            found_dest <- FALSE
                            while ( found_dest==FALSE ) {
                                dest_g <- ersample( (1:N_groups)[-g] , size=1 )
                                if ( any(is.na(pop[t,dest_g,,1])) )
                                    found_dest <- TRUE
                            }
                            dest_i <- ersample( which(is.na(pop[t,dest_g,,1])) , size=1 )
                            # remove from native group
                            pop[t,g,i,1] <- NA
                            pop[t,g,i,2] <- NA
                            # add to destination group
                            still_useful <- sample( 0:1 , size=1 , prob=c(1-v,v) )
                            pop[t,dest_g,dest_i,1] <- ff[1] * still_useful
                            pop[t,dest_g,dest_i,2] <- ff[2]
                        }#m
                    }#!is.na
                }#i
            }#g
        }# migration

    }#t

    # prep analytical quantities

    # compute average entropy (Shannon-Wiener) within groups
    # this measures within-group diversity
    H <- function(p) -sum(ifelse(p==0,0,p*log(p)))
    groupH <- function(x) {
        y <- table(x)
        p <- y/sum(y)
        H(p)
    }
    # function to compute average across group at time t
    EgH <- function(t) {
        EH <- sapply( 1:N_groups , function(g) groupH(pop[t,g,,1]) )
        mean(EH)
    }
    # compute across time
    EHwithin <- sapply( 1:tmax , function(t) EgH(t) )

    # compute average entropy across groups
    # do this by computing TOTAL population entropy and subtract EHwithin
    # Htot = EHwithin + Hbetween
    # Hbetween = Htot - EHwithin
    Htot <- sapply( 1:tmax , function(t) groupH( as.vector(pop[t,,,1]) ) )
    Hbetween <- Htot - EHwithin

    # compute frequency adaptive behavior in population
    q <- sapply( 1:tmax , function(t) mean( pop[t,,,1]>0 , na.rm=TRUE ) )

    # census size
    Ng <- array(NA,dim=c(tmax,N_groups))
    for ( t in 1:tmax )
        for ( g in 1:N_groups ) {
            Ng[t,g] <- sum( !is.na(pop[t,g,,1]) )
        }

    # compute average number of unique adaptive behaviors within groups
    ENwithin <- sapply( 1:tmax , function(t) 
        mean( 
            sapply( 1:N_groups , 
                function(g) { 
                    bb <- pop[t,g,,1]
                    bb <- bb[bb>0]
                    length(erunique(bb))
                } ) ) )

    ################
    # result
    result <- list(
            pop=pop,
            EHwithin=EHwithin,
            Hbetween=Hbetween,
            ENwithin=ENwithin,
            q=q,
            N=Ng
        )

    return(result)

}# end of sim_accumu_flex

# test code follows
if ( FALSE ) {

r <- NULL
r <- sim_accum_flex( tmax=1000 , N=100 , N_groups=10 , s=0.01 , f=0.1 , d=0.9 , e=0 , u=0 , v=0 , m=0 , death=0.05 , b=0.08 , K=0.001 , n_behaviors=1e4 , n_teachers=1 )

# library(rethinking)
# blank(h=1.8,w=2.5)
par(mfrow=c(2,3))
plot( r$q , xlab="time" , ylab="freq adaptive (population)" , ylim=c(0,1) )
plot( r$EHwithin , xlab="time" , ylab="diversity within groups" )
plot( r$Hbetween , xlab="time" , ylab="diversity between groups" )
plot( r$ENwithin , xlab="time" , ylab="number techniques per group" )
plot( apply(r$N,1,mean) , ylab="average group size" , xlab="time" )
mean_age <- sapply( 1:dim(r$pop)[1] , function(t) mean(r$pop[t,,,2],na.rm=TRUE) )
plot( mean_age , ylab="average age" , xlab="time" )

# expected group size
# log(b - death + 1)/K
log(0.08 - 0.05 + 1)/0.001

t( r$pop[800,,,1] )

}
