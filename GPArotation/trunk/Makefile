PREREQS=

include ../Makevars

default:	$(DEFAULT)

$(CODETEST): $(PREREQS:%=../%) $(PARTS) $(CODESCAN) 
	${codetestMacro}

include ../Makerules
