library( Epi )
library( popEpi )
library( survival )
data( DMlate )

dml <- Lexis( entry=list(Per=dodm, Age=dodm-dobth, DMdur=0 ),
              exit=list(Per=dox),
              exit.status=factor(!is.na(dodth),labels=c("DM","Dead")),
              data=DMlate )

dmi <- cutLexis( dml, cut = dml$doins,
                 pre = "DM",
                 new.state = "Ins",
                 new.scale = "tsI",
                 split.states = FALSE )

summary( dmi, t=T )
boxes( dmi, boxpos=TRUE, scale.R=1000, hmult=1.5 )
dms <- splitMulti( dmi, Age=seq(0,100,1/6) )

dms$tsI <- ifelse( is.na(dms$tsI), 0, dms$tsI )
nk <- 5
( a.kn <- c(40, with( subset(dms,lex.Xst=="Dead"),
                      quantile(Age+lex.dur,(1:nk-0.5)/nk) ) ) )
( i.kn <- c( 0, with( subset(dms,lex.Xst=="Dead" & lex.Cst=="Ins"),
                      quantile(tsI+lex.dur,1:(nk-1)/nk ) ) ) )

pm <- glm( cbind(lex.Xst=="Dead",lex.dur) ~ Ns(Age,knots=a.kn)
           + Ns(tsI,knots=i.kn)
           + lex.Cst + sex,
           family=poisreg, data = dms )
cm <- coxph( Surv(Age,Age+lex.dur,lex.Xst=="Dead") ~
               Ns(tsI,knots=i.kn) + factor(lex.Cst) + sex,
             data = dms )
round( cbind( ci.exp( pm, subset=c("tsI","Ins","F") ),
              + ci.exp( cm, subset=c("tsI","Ins","F") ) ), 3 )

nd <- data.frame( tsI=seq(0,15,length=151), lex.Cst="Ins", sex="M" )
nr <- data.frame( tsI= 2 , lex.Cst="Ins", sex="M" )
ppr <- ci.exp( pm, list(nd,nr), xvars="Age" )
cpr <- ci.exp( cm, list(nd,nr) )
par( mar=c(3,3,1,1), mgp=c(3,1,0)/1.6, las=1, bty="n" )
matshade( nd$tsI, cbind(ppr,cpr), lty=c(1,2), log="y", plot=T,
            ylab="Mortality RR", xlab="Time since insulin" )
abline( v=2, h=1, lty=3 )

nd <- data.frame( tsI=seq(0,15,length=151), lex.Cst="Ins", sex="M" )
nr <- data.frame( tsI= 0 , lex.Cst="DM" , sex="M" )
ppr <- ci.exp( pm, list(nd,nr), xvars="Age" )
cpr <- ci.exp( cm, list(nd,nr) )
par( mar=c(3,3,1,1), mgp=c(3,1,0)/1.6, las=1, bty="n" )
matshade( nd$tsI, cbind(ppr,cpr), lty=c(1,2), log="y", plot=T,
          ylab="Mortality RR, Ins vs. noIns", xlab="Time since insulin" )

