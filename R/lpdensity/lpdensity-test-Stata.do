clear all
set more off
import delimited "/Users/xinweima/Dropbox/0_Research/01_2017_LocPolDensity/Package_Stata_lpdensity/lpdensity_data.csv", encoding(ISO-8859-2) clear 

//--------------------------------------------------------------------------------
// Density estimation check

lpbwdensity v1, grid(g1) bwselect("mse-rot") 
lpbwdensity v2, grid(g2) bwselect("mse-rot") 
lpbwdensity v3, grid(g1) bwselect("mse-rot") 
lpbwdensity v4, grid(g2) bwselect("mse-rot") 


lpbwdensity v1, grid(g1) bwselect("mse-rot") p(1)
lpbwdensity v2, grid(g2) bwselect("mse-rot") p(1)
lpbwdensity v3, grid(g1) bwselect("mse-rot") p(1)
lpbwdensity v4, grid(g2) bwselect("mse-rot") p(1)

lpbwdensity v1, grid(g1) bwselect("mse-rot") p(3)
lpbwdensity v2, grid(g2) bwselect("mse-rot") p(3)
lpbwdensity v3, grid(g1) bwselect("mse-rot") p(3)
lpbwdensity v4, grid(g2) bwselect("mse-rot") p(3)


lpbwdensity v1, grid(g1) bwselect("imse-rot")
lpbwdensity v2, grid(g2) bwselect("imse-rot")
lpbwdensity v3, grid(g1) bwselect("imse-rot")
lpbwdensity v4, grid(g2) bwselect("imse-rot")


lpbwdensity v1, grid(g1) bwselect("imse-rot") p(1)
lpbwdensity v2, grid(g2) bwselect("imse-rot") p(1)
lpbwdensity v3, grid(g1) bwselect("imse-rot") p(1)
lpbwdensity v4, grid(g2) bwselect("imse-rot") p(1)

lpbwdensity v1, grid(g1) bwselect("imse-rot") p(3)
lpbwdensity v2, grid(g2) bwselect("imse-rot") p(3)
lpbwdensity v3, grid(g1) bwselect("imse-rot") p(3)
lpbwdensity v4, grid(g2) bwselect("imse-rot") p(3)


lpbwdensity v1, grid(g1) bwselect("mse-dpi")
lpbwdensity v2, grid(g2) bwselect("mse-dpi") 
lpbwdensity v3, grid(g1) bwselect("mse-dpi") 
lpbwdensity v4, grid(g2) bwselect("mse-dpi") 


lpbwdensity v1, grid(g1) bwselect("mse-dpi") p(1)
lpbwdensity v2, grid(g2) bwselect("mse-dpi") p(1)
lpbwdensity v3, grid(g1) bwselect("mse-dpi") p(1)
lpbwdensity v4, grid(g2) bwselect("mse-dpi") p(1)

lpbwdensity v1, grid(g1) bwselect("mse-dpi") p(3)
lpbwdensity v2, grid(g2) bwselect("mse-dpi") p(3)
lpbwdensity v3, grid(g1) bwselect("mse-dpi") p(3)
lpbwdensity v4, grid(g2) bwselect("mse-dpi") p(3)


lpbwdensity v1, grid(g1) bwselect("imse-dpi")
lpbwdensity v2, grid(g2) bwselect("imse-dpi")
lpbwdensity v3, grid(g1) bwselect("imse-dpi")
lpbwdensity v4, grid(g2) bwselect("imse-dpi")


lpbwdensity v1, grid(g1) bwselect("imse-dpi") p(1)
lpbwdensity v2, grid(g2) bwselect("imse-dpi") p(1)
lpbwdensity v3, grid(g1) bwselect("imse-dpi") p(1)
lpbwdensity v4, grid(g2) bwselect("imse-dpi") p(1)

lpbwdensity v1, grid(g1) bwselect("imse-dpi") p(3)
lpbwdensity v2, grid(g2) bwselect("imse-dpi") p(3)
lpbwdensity v3, grid(g1) bwselect("imse-dpi") p(3)
lpbwdensity v4, grid(g2) bwselect("imse-dpi") p(3)

//--------------------------------------------------------------------------------
// DF estimation check

lpbwdensity v1, grid(g1) bwselect("mse-rot") v(0) p(1)
lpbwdensity v2, grid(g2) bwselect("mse-rot") v(0) p(1)
lpbwdensity v3, grid(g1) bwselect("mse-rot") v(0) p(1)
lpbwdensity v4, grid(g2) bwselect("mse-rot") v(0) p(1)


lpbwdensity v1, grid(g1) bwselect("mse-rot") v(0) p(0)
lpbwdensity v2, grid(g2) bwselect("mse-rot") v(0) p(0)
lpbwdensity v3, grid(g1) bwselect("mse-rot") v(0) p(0)
lpbwdensity v4, grid(g2) bwselect("mse-rot") v(0) p(0)

lpbwdensity v1, grid(g1) bwselect("mse-rot") v(0) p(2)
lpbwdensity v2, grid(g2) bwselect("mse-rot") v(0) p(2)
lpbwdensity v3, grid(g1) bwselect("mse-rot") v(0) p(2)
lpbwdensity v4, grid(g2) bwselect("mse-rot") v(0) p(2)


lpbwdensity v1, grid(g1) bwselect("imse-rot") v(0) p(1)
lpbwdensity v2, grid(g2) bwselect("imse-rot") v(0) p(1)
lpbwdensity v3, grid(g1) bwselect("imse-rot") v(0) p(1)
lpbwdensity v4, grid(g2) bwselect("imse-rot") v(0) p(1)


lpbwdensity v1, grid(g1) bwselect("imse-rot") v(0) p(0)
lpbwdensity v2, grid(g2) bwselect("imse-rot") v(0) p(0)
lpbwdensity v3, grid(g1) bwselect("imse-rot") v(0) p(0)
lpbwdensity v4, grid(g2) bwselect("imse-rot") v(0) p(0)

lpbwdensity v1, grid(g1) bwselect("imse-rot") v(0) p(2)
lpbwdensity v2, grid(g2) bwselect("imse-rot") v(0) p(2)
lpbwdensity v3, grid(g1) bwselect("imse-rot") v(0) p(2)
lpbwdensity v4, grid(g2) bwselect("imse-rot") v(0) p(2)


lpbwdensity v1, grid(g1) bwselect("mse-dpi") v(0) p(1) 
lpbwdensity v2, grid(g2) bwselect("mse-dpi") v(0) p(1)
lpbwdensity v3, grid(g1) bwselect("mse-dpi") v(0) p(1) 
lpbwdensity v4, grid(g2) bwselect("mse-dpi") v(0) p(1)


lpbwdensity v1, grid(g1) bwselect("mse-dpi") v(0) p(0)
lpbwdensity v2, grid(g2) bwselect("mse-dpi") v(0) p(0)
lpbwdensity v3, grid(g1) bwselect("mse-dpi") v(0) p(0)
lpbwdensity v4, grid(g2) bwselect("mse-dpi") v(0) p(0)

lpbwdensity v1, grid(g1) bwselect("mse-dpi") v(0) p(2)
lpbwdensity v2, grid(g2) bwselect("mse-dpi") v(0) p(2)
lpbwdensity v3, grid(g1) bwselect("mse-dpi") v(0) p(2)
lpbwdensity v4, grid(g2) bwselect("mse-dpi") v(0) p(2)


lpbwdensity v1, grid(g1) bwselect("imse-dpi") v(0) p(1)
lpbwdensity v2, grid(g2) bwselect("imse-dpi") v(0) p(1)
lpbwdensity v3, grid(g1) bwselect("imse-dpi") v(0) p(1)
lpbwdensity v4, grid(g2) bwselect("imse-dpi") v(0) p(1)


lpbwdensity v1, grid(g1) bwselect("imse-dpi") v(0) p(0)
lpbwdensity v2, grid(g2) bwselect("imse-dpi") v(0) p(0)
lpbwdensity v3, grid(g1) bwselect("imse-dpi") v(0) p(0)
lpbwdensity v4, grid(g2) bwselect("imse-dpi") v(0) p(0)

lpbwdensity v1, grid(g1) bwselect("imse-dpi") v(0) p(2)
lpbwdensity v2, grid(g2) bwselect("imse-dpi") v(0) p(2)
lpbwdensity v3, grid(g1) bwselect("imse-dpi") v(0) p(2)
lpbwdensity v4, grid(g2) bwselect("imse-dpi") v(0) p(2)


//--------------------------------------------------------------------------------
// Density estimation check

lpdensity v1, grid(g1) bwselect("mse-rot") 
lpdensity v2, grid(g2) bwselect("mse-rot")
lpdensity v3, grid(g1) bwselect("mse-rot") 
lpdensity v4, grid(g2) bwselect("mse-rot")


lpdensity v1, grid(g1) bwselect("mse-rot") p(1)
lpdensity v2, grid(g2) bwselect("mse-rot") p(1)
lpdensity v3, grid(g1) bwselect("mse-rot") p(1)
lpdensity v4, grid(g2) bwselect("mse-rot") p(1)

lpdensity v1, grid(g1) bwselect("mse-rot") p(3)
lpdensity v2, grid(g2) bwselect("mse-rot") p(3)
lpdensity v3, grid(g1) bwselect("mse-rot") p(3)
lpdensity v4, grid(g2) bwselect("mse-rot") p(3)


lpdensity v1, grid(g1) bwselect("imse-rot")
lpdensity v2, grid(g2) bwselect("imse-rot")
lpdensity v3, grid(g1) bwselect("imse-rot")
lpdensity v4, grid(g2) bwselect("imse-rot")


lpdensity v1, grid(g1) bwselect("imse-rot") p(1)
lpdensity v2, grid(g2) bwselect("imse-rot") p(1)
lpdensity v3, grid(g1) bwselect("imse-rot") p(1)
lpdensity v4, grid(g2) bwselect("imse-rot") p(1)

lpdensity v1, grid(g1) bwselect("imse-rot") p(3)
lpdensity v2, grid(g2) bwselect("imse-rot") p(3)
lpdensity v3, grid(g1) bwselect("imse-rot") p(3)
lpdensity v4, grid(g2) bwselect("imse-rot") p(3)


lpdensity v1, grid(g1) bwselect("mse-dpi")
lpdensity v2, grid(g2) bwselect("mse-dpi") 
lpdensity v3, grid(g1) bwselect("mse-dpi")
lpdensity v4, grid(g2) bwselect("mse-dpi")


lpdensity v1, grid(g1) bwselect("mse-dpi") p(1)
lpdensity v2, grid(g2) bwselect("mse-dpi") p(1)
lpdensity v3, grid(g1) bwselect("mse-dpi") p(1)
lpdensity v4, grid(g2) bwselect("mse-dpi") p(1)

lpdensity v1, grid(g1) bwselect("mse-dpi") p(3)
lpdensity v2, grid(g2) bwselect("mse-dpi") p(3)
lpdensity v3, grid(g1) bwselect("mse-dpi") p(3)
lpdensity v4, grid(g2) bwselect("mse-dpi") p(3)


lpdensity v1, grid(g1) bwselect("imse-dpi")
lpdensity v2, grid(g2) bwselect("imse-dpi")
lpdensity v3, grid(g1) bwselect("imse-dpi")
lpdensity v4, grid(g2) bwselect("imse-dpi")


lpdensity v1, grid(g1) bwselect("imse-dpi") p(1)
lpdensity v2, grid(g2) bwselect("imse-dpi") p(1)
lpdensity v3, grid(g1) bwselect("imse-dpi") p(1)
lpdensity v4, grid(g2) bwselect("imse-dpi") p(1)

lpdensity v1, grid(g1) bwselect("imse-dpi") p(3)
lpdensity v2, grid(g2) bwselect("imse-dpi") p(3)
lpdensity v3, grid(g1) bwselect("imse-dpi") p(3)
lpdensity v4, grid(g2) bwselect("imse-dpi") p(3)

//--------------------------------------------------------------------------------
// DF estimation check

lpdensity v1, grid(g1) bwselect("mse-rot") v(0) p(1)
lpdensity v2, grid(g2) bwselect("mse-rot") v(0) p(1)
lpdensity v3, grid(g1) bwselect("mse-rot") v(0) p(1)
lpdensity v4, grid(g2) bwselect("mse-rot") v(0) p(1)


lpdensity v1, grid(g1) bwselect("mse-rot") v(0) p(0)
lpdensity v2, grid(g2) bwselect("mse-rot") v(0) p(0)
lpdensity v3, grid(g1) bwselect("mse-rot") v(0) p(0)
lpdensity v4, grid(g2) bwselect("mse-rot") v(0) p(0)

lpdensity v1, grid(g1) bwselect("mse-rot") v(0) p(2)
lpdensity v2, grid(g2) bwselect("mse-rot") v(0) p(2)
lpdensity v3, grid(g1) bwselect("mse-rot") v(0) p(2)
lpdensity v4, grid(g2) bwselect("mse-rot") v(0) p(2)


lpdensity v1, grid(g1) bwselect("imse-rot") v(0) p(1)
lpdensity v2, grid(g2) bwselect("imse-rot") v(0) p(1)
lpdensity v3, grid(g1) bwselect("imse-rot") v(0) p(1)
lpdensity v4, grid(g2) bwselect("imse-rot") v(0) p(1)


lpdensity v1, grid(g1) bwselect("imse-rot") v(0) p(0)
lpdensity v2, grid(g2) bwselect("imse-rot") v(0) p(0)
lpdensity v3, grid(g1) bwselect("imse-rot") v(0) p(0)
lpdensity v4, grid(g2) bwselect("imse-rot") v(0) p(0)

lpdensity v1, grid(g1) bwselect("imse-rot") v(0) p(2)
lpdensity v2, grid(g2) bwselect("imse-rot") v(0) p(2)
lpdensity v3, grid(g1) bwselect("imse-rot") v(0) p(2)
lpdensity v4, grid(g2) bwselect("imse-rot") v(0) p(2)


lpdensity v1, grid(g1) bwselect("mse-dpi") v(0) p(1) 
lpdensity v2, grid(g2) bwselect("mse-dpi") v(0) p(1)
lpdensity v3, grid(g1) bwselect("mse-dpi") v(0) p(1) 
lpdensity v4, grid(g2) bwselect("mse-dpi") v(0) p(1)

lpdensity v1, grid(g1) bwselect("mse-dpi") v(0) p(0)
lpdensity v2, grid(g2) bwselect("mse-dpi") v(0) p(0)
lpdensity v3, grid(g1) bwselect("mse-dpi") v(0) p(0)
lpdensity v4, grid(g2) bwselect("mse-dpi") v(0) p(0)

lpdensity v1, grid(g1) bwselect("mse-dpi") v(0) p(2)
lpdensity v2, grid(g2) bwselect("mse-dpi") v(0) p(2)
lpdensity v3, grid(g1) bwselect("mse-dpi") v(0) p(2)
lpdensity v4, grid(g2) bwselect("mse-dpi") v(0) p(2)


lpdensity v1, grid(g1) bwselect("imse-dpi") v(0) p(1)
lpdensity v2, grid(g2) bwselect("imse-dpi") v(0) p(1)
lpdensity v3, grid(g1) bwselect("imse-dpi") v(0) p(1)
lpdensity v4, grid(g2) bwselect("imse-dpi") v(0) p(1)


lpdensity v1, grid(g1) bwselect("imse-dpi") v(0) p(0)
lpdensity v2, grid(g2) bwselect("imse-dpi") v(0) p(0)
lpdensity v3, grid(g1) bwselect("imse-dpi") v(0) p(0)
lpdensity v4, grid(g2) bwselect("imse-dpi") v(0) p(0)

lpdensity v1, grid(g1) bwselect("imse-dpi") v(0) p(2)
lpdensity v2, grid(g2) bwselect("imse-dpi") v(0) p(2)
lpdensity v3, grid(g1) bwselect("imse-dpi") v(0) p(2)
lpdensity v4, grid(g2) bwselect("imse-dpi") v(0) p(2)
