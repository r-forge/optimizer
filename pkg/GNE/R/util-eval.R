evalwitharglist <- function(f, x, addarg)
{
	if(length(addarg) == 0)
	{	
		f(x)
	}else
	{	
		do.call(f, c(list(x), addarg) )
	}	
}
