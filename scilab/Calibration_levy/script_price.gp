set title 'Price'
set xlabel 'Strike/Spot'
set ylabel 'Maturity'
splot "./Data/MarketData_withvol.dat" title 'market price', "./Data/ModelData.dat" title 'fit price'
set terminal postscript eps color
set output "price.eps"
replot
