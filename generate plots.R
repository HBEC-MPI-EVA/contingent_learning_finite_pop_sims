# function to generate plots for paper
# can be set to run sims -or- use prepared sim data from function

make_accum_dat <- function( 
    sims=NULL , 
    N_groups = 16,
    s = 0.01,
    d = 1,
    e = 0.1,
    u = 0.01,
    v = 1,
    m = 0.01,
    death = 0.05,
    b = 0.1,
    K = log(0.1-0.05+1)/20,
    f = 1, 
    n_behaviors = 1e4,
    n_teachers = 1,
    tmax = 2000,
    end_range = 500,
    n_sims = 1,
    cores = 3,
    preschedule = FALSE ,
    ... ) 
{

    the_call <- match.call()
    mnames <- names(the_call)[-1]
    # list of arguments to recognize as graphed parameters
    white_list <- c( "s","d","e","u","v","m","death","b","K","f","n_behaviors","n_teachers" )
    mnames <- mnames[ mnames %in% white_list ]
    print(mnames)

    library(parallel)
    library(rethinking)

    # set up list of parameter values
    seq <- expand.grid(
        N_groups = N_groups,
        s = s,
        d = d,
        e = e,
        u = u,
        v = v,
        m = m,
        f = f,
        death = death,
        b = b,
        K = K, # gives different average group sizes
        n_behaviors = n_behaviors,
        n_teachers = n_teachers
    )

    # set up function wrapper that calls simulation and cleans results
    library(compiler)
    sim_cmp <- cmpfun(sim_accum_flex)

    t2 <- tmax-end_range
    ff <- function(i) {
        p <- seq[i,]
        #print(p)
        ss <- with( p , sim_cmp( 
            tmax=tmax , 
            N=100 ,
            N_init=log(b-death+1)/K , 
            N_groups=N_groups , 
            s=s , 
            d=d , 
            e=e , 
            u=u , 
            v=v , 
            m=m , 
            f=p$f ,
            death=death , b=b , K=K , 
            n_behaviors=n_behaviors , n_teachers=n_teachers ) )
        r <- list(
            EHwithin = mean(ss$EHwithin[t2:tmax]),
            Hbetween = mean(ss$Hbetween[t2:tmax]),
            ENwithin = mean(ss$ENwithin[t2:tmax]),
            q = mean(ss$q[t2:tmax]),
            N = mean(ss$N[t2:tmax])
        )
        return(r)
    }
    f_cmp <- cmpfun(ff)

    # dispatch to cores
    num_cores <- cores
    if (n_sims>1) {
        # more than one replicate for each combination of parameters
        # so duplicate the rows with rbind
        seq_orig <- seq
        for ( i in 1:n_sims ) seq <- rbind(seq,seq_orig)
    }
    nn <- nrow(seq)
    result <- mclapply( 
        1:nrow(seq) ,
        function(i) {
            print(concat(i,"/",nn))
            return(f_cmp(i))
        } , 
        mc.cores=num_cores , mc.preschedule=preschedule )

    # simplify from list to array
    result2 <- t( simplify2array(result) )
    result3 <- as.data.frame(matrix(unlist(result), nrow=length(unlist(result[1])), byrow=FALSE))
    result3 <- t(result3)
    colnames(result3) <- colnames(result2)
    rownames(result3) <- NULL

    # merge with parameter values
    result4 <- cbind( seq[1:nrow(result3),] , result3 )
    attr(result4,"parnames") <- mnames
    attr(result4,"call") <- the_call

    invisible(result4)

}#

make_accum_plot <- function( sim_set , y="q" , x=NULL , z=NULL , ylim=NULL , xlim=NULL , new=TRUE , logx=FALSE , idpts=TRUE , xlab , ylab , cex=0.7 , na.rm=TRUE , ... ) {
    pars <- attr(sim_set,"parnames")
    if ( is.null(x) ) {
        x <- pars[1]
        z <- pars[2]
    }
    if ( y=="H" ) {
        sim_set$H <- sim_set$EHwithin + sim_set$Hbetween
    }
    if ( y=="P" ) {
        sim_set$P <- sim_set$Hbetween / ( sim_set$EHwithin + sim_set$Hbetween )
    }
    if ( is.null(ylim) ) ylim <- range(sim_set[[y]])
    if ( new==TRUE ) blank()
    if ( is.null(xlim) ) xlim <- range(sim_set[[x]])
    if ( logx==TRUE ) xlim <- log(xlim)
    if ( missing(xlab) ) xlab <- x
    if ( missing(ylab) ) ylab <- y
    plot( NULL , xlim=xlim , ylim=ylim , xlab=xlab , ylab=ylab , ... )
    xu <- unique(sim_set[[x]])
    zu <- unique(sim_set[[z]])
    #axis( 1 , at=xu , labels=xu )
    xp <- xu
    if ( logx==TRUE ) xp <- log(xu)
    for ( zi in zu ) {
        #yu <- sim_set[[y]][sim_set[[z]]==zi]
        yu <- xp # just need dimensions
        for ( i in 1:length(xp) ) {
            yu[i] <- mean( sim_set[[y]][ sim_set[[z]]==zi & sim_set[[x]]==xp[i] ] , na.rm=na.rm )
        }
        lines( xp , yu , ... )
        points( xp , yu , pch=16 , ... )
    }#zi
    if ( idpts==TRUE ) {
        labs <- paste( z , "=" , sim_set[[z]] )
        identify( sim_set[[x]] , sim_set[[y]] , labels=labs , atpen=TRUE , cex=cex )
    }
}

# example use

zz <- make_accum_dat( 
    m=c(0,0.01,0.04,0.1) , u=c(0,0.01,0.05,0.1) , 
    s=0.001 , K = log(0.1-0.05+1)/20 , b=0.1 , death=0.05 , n_teachers=1 , 
    tmax=5000 , end_range=1000 )

make_accum_plot( zz , y="q" , x="m" , z="u" , ylim=c(0,1) , logx=FALSE , idpts=FALSE )

############################ impact of migration
# q * m * u [and H,p]

zz <- make_accum_dat( 
    m=c(0,0.01,0.04,0.1) , u=c(0.001,0.01,0.1) , 
    s=0.001 , K = log(0.1-0.01+1)/10 , b=0.1 , death=0.01 , n_teachers=1 , 
    tmax=5000 , end_range=2000 , cores=60 , n_sims=10 )
save( zz , file="zz_qmu.rda" )

# load( "zz_qmu.rda" )

blank(w=3,ex=0.7)
par(mfrow=c(1,3))
make_accum_plot( zz , y="q" , x="m" , z="u" , ylim=c(0,1) , logx=FALSE , idpts=TRUE , new=FALSE , ylab="proportion successful" )
make_accum_plot( zz , y="H" , x="m" , z="u" , ylim=c(0,2.5) , logx=FALSE , idpts=TRUE , new=FALSE , ylab="total diversity (H)" )
make_accum_plot( zz , y="P" , x="m" , z="u" , ylim=c(0,1) , logx=FALSE , idpts=TRUE , new=FALSE , ylab="proportion diversity between groups" )

############################ impact of group size AND migration
# q * N * m [and H,p]

zz <- make_accum_dat( 
    m=c(0.001,0.01,0.1) , K = log(0.1-0.01+1)/c(10,20,50) , 
    u=0.01 , s=0.001 , b=0.1 , death=0.01 , n_teachers=1 , 
    tmax=5000 , end_range=2000 , cores=60 , n_sims=10 )
save( zz , file="zz_qNm.rda" )

# load( "zz_qNm.rda" )

blank(w=3,ex=0.7)
par(mfrow=c(1,3))
make_accum_plot( zz , y="q" , x="K" , z="m" , ylim=c(0,1) , logx=FALSE , idpts=TRUE , new=FALSE , ylab="proportion successful" , xaxt="n" , xlab="N" )
axis( 1 , at=log(0.1-0.01+1)/c(10,20,50) , labels=c(10,20,50) )
make_accum_plot( zz , y="H" , x="K" , z="m" , ylim=c(0,3) , logx=FALSE , idpts=TRUE , new=FALSE , ylab="total diversity (H)" , xaxt="n" , xlab="N" )
axis( 1 , at=log(0.1-0.01+1)/c(10,20,50) , labels=c(10,20,50) )
make_accum_plot( zz , y="P" , x="K" , z="m" , ylim=c(0,1) , logx=FALSE , idpts=TRUE , new=FALSE , ylab="proportion diversity between groups" , xaxt="n" , xlab="N" )
axis( 1 , at=log(0.1-0.01+1)/c(10,20,50) , labels=c(10,20,50) )

############################ impact of number of teachers
# q * N * m with 10 teachers [and H,p]

zz <- make_accum_dat( 
    m=c(0.001,0.01,0.1) , K = log(0.1-0.01+1)/c(10,20,50) , 
    u=0.01 , s=0.001 , b=0.1 , death=0.01 , n_teachers=10 , 
    tmax=5000 , end_range=2000 , cores=60 , n_sims=10 )
save( zz , file="zz_qNm10t.rda" )

# load( "zz_qNm10t.rda" )

blank(w=3,ex=0.7)
par(mfrow=c(1,3))
make_accum_plot( zz , y="q" , x="K" , z="m" , ylim=c(0,1) , logx=FALSE , idpts=TRUE , new=FALSE , ylab="proportion successful" , xaxt="n" , xlab="N" )
axis( 1 , at=log(0.1-0.01+1)/c(10,20,50) , labels=c(10,20,50) )
make_accum_plot( zz , y="H" , x="K" , z="m" , ylim=c(0,3) , logx=FALSE , idpts=TRUE , new=FALSE , ylab="total diversity (H)" , xaxt="n" , xlab="N" )
axis( 1 , at=log(0.1-0.01+1)/c(10,20,50) , labels=c(10,20,50) )
make_accum_plot( zz , y="P" , x="K" , z="m" , ylim=c(0,1) , logx=FALSE , idpts=TRUE , new=FALSE , ylab="proportion diversity between groups" , xaxt="n" , xlab="N" )
axis( 1 , at=log(0.1-0.01+1)/c(10,20,50) , labels=c(10,20,50) )

############################
# q * s * m [and H,p]

zz <- make_accum_dat( 
    s=c(0.001,0.005,0.01,0.05,0.1) , m=c(0.001,0.01,0.1) ,
    u=0.01 , K = log(0.1-0.01+1)/10 , b=0.1 , death=0.01 , n_teachers=1 , 
    tmax=5000 , end_range=2000 , cores=60 , n_sims=10 )
save( zz , file="zz_qsm.rda" )

# load( "zz_qsm.rda" )

blank(w=3,ex=0.7)
par(mfrow=c(1,3))
make_accum_plot( zz , y="q" , x="s" , z="m" , ylim=c(0,1) , logx=FALSE , idpts=TRUE , new=FALSE , ylab="proportion successful" , xlab="innovation rate (s)" )
make_accum_plot( zz , y="H" , x="s" , z="m" , ylim=c(0,3) , logx=FALSE , idpts=TRUE , new=FALSE , ylab="total diversity (H)" , xlab="innovation rate (s)" )
make_accum_plot( zz , y="P" , x="s" , z="m" , ylim=c(0,1) , logx=FALSE , idpts=TRUE , new=FALSE , ylab="proportion diversity between groups" , xlab="innovation rate (s)" )

############################
# q * s * e [and H,p]

zz <- make_accum_dat( 
    s=c(0.001,0.01,0.1) , e=c(0.1,0.2,0.5) , m=0.01 ,
    u=0.01 , K = log(0.1-0.01+1)/10 , b=0.1 , death=0.01 , n_teachers=1 , 
    tmax=5000 , end_range=2000 , cores=60 , n_sims=10 )
save( zz , file="zz_qse.rda" )

# load( "zz_qse.rda" )

blank(w=3,ex=0.7)
par(mfrow=c(1,3))
make_accum_plot( zz , y="q" , x="s" , z="e" , ylim=c(0,1) , logx=FALSE , idpts=TRUE , new=FALSE , ylab="proportion successful" , xlab="innovation rate (s)" )
make_accum_plot( zz , y="H" , x="s" , z="e" , ylim=c(0,3) , logx=FALSE , idpts=TRUE , new=FALSE , ylab="total diversity (H)" , xlab="innovation rate (s)" )
make_accum_plot( zz , y="P" , x="s" , z="e" , ylim=c(0,1) , logx=FALSE , idpts=TRUE , new=FALSE , ylab="proportion diversity between groups" , xlab="innovation rate (s)" )

