ilib_name="libmanu"

files=["implied.o";"implied_volatility.o"];

libs=[];

table=["implied_volatility","implied"];

ilib_build(ilib_name,table,files,libs);
