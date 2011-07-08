# The built-in widget "gtext" in gWidgets is just too slow,
# so here we define a new one named "gtextbox", which proves to be faster.
setClass("myTextBox", representation(textbox = "gGroupRGtk",
				textbuf = "GtkTextBuffer",
				textview = "GtkTextView"));
# The constructor of "gtextbox"
gtextbox = function(text = "", wrap = TRUE, font.attr = NULL,
		            container = NULL, margin = c(10, 10),
					frame = FALSE, frameLabel = "Code", ...)
{
	hBox = gtkHBoxNew();
	alignBox = gtkAlignmentNew(0.5, 0, 1, 1);
	alignBox$setPadding(0, 0, 20, 20);
	textBuf = gtkTextBufferNew();
	textView = gtkTextView(buffer = textBuf, show = TRUE);
	if(wrap) textView$SetWrapMode(GtkWrapMode["word"]);
	textView$setLeftMargin(margin[1]);
	textView$setRightMargin(margin[2]);
	hBox$packStart(alignBox);
	if(frame)
	{
		frameBox = gtkFrameNew(frameLabel);
		alignBox$add(frameBox);
		frameBox$add(textView);
	}else{
		alignBox$add(textView);
	}
	textBox = as.gWidgetsRGtk2(hBox);
	obj = new("myTextBox", textbox = textBox, textview = textView,
           textbuf = textBuf);
	
	if(is.null(font.attr)) font.attr = c(family = "sans", style = "normal",
				                         weight = "normal", size = 11,
				                         color = "black");
	setfont.myTextBox(obj, value = font.attr);
	setvalue.myTextBox(obj, text);
	if(!is.null(container)) add(container, textBox);
	return(obj);
}
# Member functions
getvalue.myTextBox = function(obj)
{
	buf = obj@textbuf;
	startIter = buf$getStartIter()$iter;
	endIter = buf$getEndIter()$iter;
	return(buf$getText(startIter, endIter));
}
setvalue.myTextBox = function(obj, value)
{
	obj@textbuf$setText(as.character(value));
	invisible(obj);
}
getfont.myTextBox = function(obj)
{
	pangoContext = obj@textview$getPangoContext();
	fontDescrip = pangoContext$getFontDescription();
	font.attr = c(family = fontDescrip$getFamily(),
			      style = names(PangoStyle[PangoStyle == as.integer(fontDescrip$getStyle())]),
			      weight = names(PangoWeight[PangoWeight == as.integer(fontDescrip$getWeight())]),
				  size = fontDescrip$getSize() / PANGO_SCALE);
	return(font.attr);
}
setfont.myTextBox = function(obj, value)
{
	pangoContext = obj@textview$getPangoContext();
	fontDescrip = pangoContext$getFontDescription();
	if("family" %in% names(value))
	{
		fontDescrip$setFamily(value["family"]);
	}
	if("style" %in% names(value))
	{
		fontDescrip$setStyle(PangoStyle[value["style"]]);
	}
	if("weight" %in% names(value))
	{
		fontDescrip$setWeight(PangoWeight[value["weight"]]);
	}
	if("size" %in% names(value))
	{
		fontDescrip$setSize(as.integer(value["size"]) * PANGO_SCALE);
	}
	if("color" %in% names(value))
	{
		obj@textview$modifyText(GtkStateType["normal"],
				                 gdkColorParse(value["color"])$color);
	}
	obj@textview$modifyFont(fontDescrip);
	invisible(obj);
}
setbasecolor.myTextBox = function(obj, r, g, b)
{
	color = rgb(r / 65536, g / 65536, b / 65536);
	obj@textview$modifyBase(GtkStateType["normal"],
			                 gdkColorParse(color)$color);
}
setvisible.myTextBox = function(obj, value)
{
    visible(obj@textbox) = as.logical(value);
    invisible(obj);
}
# Set methods
setMethod("svalue", "myTextBox", getvalue.myTextBox);
setMethod("font", "myTextBox", getfont.myTextBox);
setReplaceMethod("svalue", "myTextBox", setvalue.myTextBox);
setReplaceMethod("font", "myTextBox", setfont.myTextBox);
setReplaceMethod("visible", "myTextBox", setvisible.myTextBox);
# End of definition of "gtextbox"


# Another widget is a label that could be "clickable"
# Link label
glinklabel = function(text = "", ...)
{
	label = gtkLabelNew();
	label$setAlignment(0, 0.5);
	label$setMarkup(sprintf("<a href=''>%s</a>", text));
	gSignalConnect(label, "activate-link", function(widget, uri, param) return(TRUE),
			NULL);
	return(as.gWidgetsRGtk2(label));
}
# End of definition of "glinklabel"