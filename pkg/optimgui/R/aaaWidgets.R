# The built-in widget "gtext" in "gWidgets" is not enough,
# so here we define a new one named "gtextbox".

setClass("gTextBox", representation(textbox = "gGroupRGtk",
                                    textbuf = "GtkTextBuffer",
                                    textview = "GtkTextView"));
# The "constructor" of "gtextbox"
gtextbox = function(text = "", wrap = TRUE, font.attr = NULL,
		            container = NULL, margin = c(10, 10),
					frame = FALSE, frameLabel = "Code", ...)
{
	hBox = gtkHBoxNew();
	alignBox = gtkAlignmentNew(0.5, 0, 1, 1);
	alignBox$setPadding(0, 0, 20, 20);
scrollWindow=gtkScrolledWindowNew();
scrollWindow$setPolicy(GtkPolicyType["automatic"], GtkPolicyType["automatic"]);
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
	obj = new("gTextBox", textbox = textBox, textview = textView,
               textbuf = textBuf);
	
	if(is.null(font.attr)) font.attr = c(family = "sans", style = "normal",
				                         weight = "normal", size = 11,
				                         color = "black");
	setfont.gTextBox(obj, value = font.attr);
	setvalue.gTextBox(obj, text);
	if(!is.null(container)) add(container, textBox);
	return(obj);
}

# Member functions
getvalue.gTextBox = function(obj)
{
	buf = obj@textbuf;
	startIter = buf$getStartIter()$iter;
	endIter = buf$getEndIter()$iter;
	return(buf$getText(startIter, endIter));
}
setvalue.gTextBox = function(obj, value)
{
	obj@textbuf$setText(as.character(value));
	invisible(obj);
}
getfont.gTextBox = function(obj)
{
	pangoContext = obj@textview$getPangoContext();
	fontDescrip = pangoContext$getFontDescription();
	font.attr = c(family = fontDescrip$getFamily(),
			      style = names(PangoStyle[PangoStyle == as.integer(fontDescrip$getStyle())]),
			      weight = names(PangoWeight[PangoWeight == as.integer(fontDescrip$getWeight())]),
				  size = fontDescrip$getSize() / PANGO_SCALE);
	return(font.attr);
}
setfont.gTextBox = function(obj, value)
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
getvisible.gTextBox = function(obj)
{
    invisible(obj@textbox@widget$getVisible());
}
setvisible.gTextBox = function(obj, value)
{
    visible(obj@textbox) = as.logical(value);
    invisible(obj);
}

# Register methods
setMethod("svalue", "gTextBox", getvalue.gTextBox);
setMethod("font", "gTextBox", getfont.gTextBox);
setMethod("visible", "gTextBox", getvisible.gTextBox);
setReplaceMethod("svalue", "gTextBox", setvalue.gTextBox);
setReplaceMethod("font", "gTextBox", setfont.gTextBox);
setReplaceMethod("visible", "gTextBox", setvisible.gTextBox);
# End of definition of "gtextbox"



# Another widget is a label that could be "clickable".
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
