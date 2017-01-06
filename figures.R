# figures

# show contour of qbar

u_list <- seq(from=0,to=0.2,length.out=20)
s_list <- seq(from=0.001,to=0.5,length.out=40)
d <- 1
e <- 0.1
q <- outer( s_list , u_list , function(s,u) d*s/( (1-u)*( d*s + e*(1-s) )+u ) )
contour( s_list , u_list , q , xlab="innovation rate (s)" , ylab="environmental change (u)" )
mtext( "1 teacher" )

# sample n teachers
# no explicit analytical solution, so solve numerically in each case
# just use optimize() to find zero point
# optimize(  )

fDq <- function(x,s,u) {
    Qt <- 1 - ( 1 - x*(1-u) )^n
    Dq <- Qt*(1-e+e*s) + (1-Qt)*d*s - x
    return(Dq)
}

u_list <- seq(from=0,to=0.2,length.out=20)
s_list <- seq(from=0.001,to=0.5,length.out=40)
d <- 1
e <- 0.1
n <- 2

qsols <- outer( s_list , u_list , function(a,b) a*b*0 )
for ( i in 1:length(s_list) )
    for ( j in 1:length(u_list) ) {
        qsols[i,j] <- uniroot( fDq , c(0,1) , s=s_list[i] , u=u_list[j] )[[1]]
    }

contour( s_list , u_list , qsols , xlab="innovation rate (s)" , ylab="environmental change (u)" )
mtext( "2 teachers" )

# overlapping generations

fDq2 <- function(x,s,u) {
    Qt <- 1 - ( 1 - x*(1-u) )^n
    Dq <- (1 - mu)*(Qt + (1 - Qt)*(Qt*(1 - e + e*s) + (1 - Qt)*d*s)) + mu*(Qt*(1 - e + e*s) + (1 - Qt)*d*s) - x
    return(Dq)
}
# uniroot(fDq2,c(0,1),s=0.01,u=0.1)

u_list <- seq(from=0,to=0.2,length.out=20)
s_list <- seq(from=0.0001,to=0.5,length.out=40)
d <- 1
e <- 0.1
n <- 1
mu <- 0.01

qsols <- outer( s_list , u_list , function(a,b) a*b*0 )
for ( i in 1:length(s_list) )
    for ( j in 1:length(u_list) ) {
        qsols[i,j] <- uniroot( fDq2 , c(0,1) , s=s_list[i] , u=u_list[j] )[[1]]
    }

contour( s_list , u_list , qsols , xlab="innovation rate (s)" , ylab="environmental change (u)" )
mtext( "overlapping generations" )

# show s by e

e_list <- seq(from=0,to=0.5,length.out=20)
s_list <- seq(from=0.0001,to=0.5,length.out=40)
d <- 1
rm(e)
u <- 0.1
n <- 1
mu <- 0.01

fDq3 <- function(x,s,e) {
    Qt <- 1 - ( 1 - x*(1-u) )^n
    Dq <- (1 - mu)*(Qt + (1 - Qt)*(Qt*(1 - e + e*s) + (1 - Qt)*d*s)) + mu*(Qt*(1 - e + e*s) + (1 - Qt)*d*s) - x
    return(Dq)
}
qsols <- outer( s_list , e_list , function(a,b) a*b*0 )
for ( i in 1:length(s_list) )
    for ( j in 1:length(u_list) ) {
        qsols[i,j] <- uniroot( fDq3 , c(0,1) , s=s_list[i] , e=e_list[j] )[[1]]
    }

contour( s_list , e_list , qsols , xlab="innovation rate (s)" , ylab="transmission error (e)" )
mtext( "overlapping generations" )
