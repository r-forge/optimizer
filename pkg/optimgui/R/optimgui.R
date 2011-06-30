optimgui = function()
{
options("guiToolkit" = "RGtk2");

# Rop file parser
"%+%" = function(x, y) sprintf("%s%s", x, y);
parseRop = function(filePath)
{
	RopFile = readLines(filePath);
	RopFile = RopFile[RopFile != ""];
	index.XML = grep("^#@", RopFile);
	index.rcode = (1:length(RopFile))[-index.XML];
	RopFile[index.XML] = gsub("^#@ ", "", RopFile[index.XML]);
	RopFile[index.rcode] = gsub("<", "&lt;", RopFile[index.rcode]);
	RopFile[index.rcode] = gsub(">", "&gt;", RopFile[index.rcode]);
	RopTree = xmlTreeParse(RopFile, asText = TRUE);
	return(RopTree);
}
writeRop = function(filePath)
{
	Rop = xmlTree("Roptimgui");
	for(i in 1:length(tabsList))
	{
		Rop$addNode("tab", attrs = c(tabname = names(tabsList)[i]),
				    close = FALSE);
		for(j in 1:length(tabsList[[i]]))
		{
            tag = names(tabsList[[i]])[j];
            textValue = svalue(tabsList[[i]][[j]]);
            if(tag == "label")
            {
                textValue = strwrap(textValue, 70);
                textValue = "\n" %+% paste("    ", textValue, sep = "",
                                           collapse = "\n") %+% "\n    ";
            }else if(tag == "rcode")
            {
                textValue = gsub("\n", "\n##Rcode##", textValue);
                textValue = "\n##Rcode##" %+% textValue %+% "\n    ";
            }
			Rop$addNode(tag, textValue);
		}
		Rop$closeTag();
	}
    z = textConnection("buffer", open = "w");
    sink(z);
	cat(saveXML(Rop$value()));
    sink();
    close(z);
    buffer = paste("#@ ", buffer, sep = "");
    buffer = gsub("#@ ##Rcode##", "", buffer);
	buffer = gsub("&lt;", "<", buffer);
	buffer = gsub("&gt;", ">", buffer);
    writeLines(buffer, filePath);
}
# Add the content of a tab to a notebook widget
addTabs = function(tabNode, noteBookWidget)
{
	tabBox = ggroup(horizontal = FALSE, container = noteBookWidget,
			        label = xmlAttrs(tabNode)["tabname"], spacing = 15);
	addWidgets = function(widgetNode, boxWidget)
	{
		if(xmlName(widgetNode) == "title")
		{
			titleLabel = glabel(xmlValue(widgetNode), container = boxWidget);
			font(titleLabel) = c(size = "xx-large");
			return(titleLabel);
		}else if(xmlName(widgetNode) == "label"){
			labelContent = gsub("\n +", " ", xmlValue(widgetNode));
			labelContent = paste(labelContent, collapse = "");
			docLabel = gtextbox(labelContent, container = boxWidget,
					            font.attr = c(family = "sans", size = 11));
			return(docLabel);
		}else if(xmlName(widgetNode) == "rcode"){
			codeText = gtextbox(xmlValue(widgetNode), container = boxWidget,
					            font.attr = c(family = "monospace", size = 11),
								frame = TRUE);
			return(codeText);
		}
	}
	return(xmlApply(tabNode, addWidgets, boxWidget = tabBox));
}
# Open Rop file
currentFile = NULL;
openRopFile = function(filePath)
{
	visible(buttonGroup) = FALSE;
	visible(noteBook) = FALSE;
	while(dispose(noteBook)){}
	path = filePath;
	Encoding(path) = "UTF-8";
	RopTree = parseRop(path);
	currentFile <<- path;
	tabsList <<- xmlApply(xmlRoot(RopTree), addTabs, noteBookWidget = noteBook);
	names(tabsList) <<- xmlSApply(xmlRoot(RopTree), xmlAttrs);
	visible(noteBook) = TRUE;
	visible(buttonGroup) = TRUE;
    enabled(mSave) = TRUE;
	return(NULL);
}


# Main Window
mainWin = gwindow("optimgui", visible = FALSE);
size(mainWin) = c(800, 600);
# Horizontal layout box
hGroup = ggroup(container = mainWin, expand = TRUE);
# Notebook widget
noteBook = gnotebook(closebuttons = FALSE, container = hGroup, expand = TRUE);
noteBook@widget@widget$modifyBg(GtkStateType["normal"],
		                         gdkColorParse("white")$color);

# Welcome page
welcomePageNew = function()
{
	welcomePage = ggroup(horizontal = FALSE, expand = TRUE);
	welcomeLabel = glabel("Create a new project or open an existing Rop file.",
			container = welcomePage);
	font(welcomeLabel) = c(size = "x-large");
	tmpBox1 = ggroup(container = welcomePage);
	tmpBox2 = ggroup(container = welcomePage);
	tmpBox3 = ggroup(container = welcomePage);
	gimage(system.file("resources", "images",  "wizard.png", package = "optimgui"), container = tmpBox1);
	wizardLabel = glinklabel("Create a new project by wizard");
	add(tmpBox1, wizardLabel);
	gimage(system.file("resources", "images",  "template.png", package = "optimgui"), container = tmpBox2);
	newLabel = glinklabel("Create a new project with template");
	add(tmpBox2, newLabel);
	gimage(system.file("resources", "images",  "open.png", package = "optimgui"), container = tmpBox3);
	openLabel = glinklabel("Open an existing Rop file");
	add(tmpBox3, openLabel);
	rtvalue = list(welcomePage = welcomePage, wizardLabel = wizardLabel,
			newLabel = newLabel, openLabel = openLabel);
	welcomePageAddEvent(rtvalue);
	return(rtvalue);
}

# Catalog page
catalogPageNew = function()
{
	catalogPage = ggroup(horizontal = FALSE, expand = TRUE);
	# FIXME
	Ropfiles = data.frame(FileName = "shobbs.Rop",
			Description = "Scaled Hobbs weeds problem");
	RopTable = gtable(Ropfiles, container = catalogPage, expand = TRUE);
	gseparator(container = catalogPage);
	tmpBox4 = ggroup(container = catalogPage);
	openRopButton = gbutton("OK", container = tmpBox4);
	rtvalue = list(catalogPage = catalogPage, RopTable = RopTable,
			openRopButton = openRopButton);
	catalogPageAddEvent(rtvalue);
	return(rtvalue);
}

# Wizard page. This part is programmed using RGtk2.
wizardPage = gtkAssistantNew(show = FALSE);
wizardPage$setTitle("Wizard");
wizardPage$setDefaultSize(640, 480);
# Welcome page
wizardWelcomePage = gtkLabelNew("In this wizard you will be asked several questions about your optimization problem,
				and in the final step a suggested template will be provided.");
wizardWelcomePage$modifyFont(pangoFontDescriptionFromString("sans 11"));
wizardWelcomePage$setAlignment(0.2, 0.1);
wizardPage$appendPage(wizardWelcomePage);
wizardPage$setPageTitle(wizardWelcomePage, "Welcome");
wizardPage$setPageComplete(wizardWelcomePage, TRUE);
# Question 1, just a demo
wizardQ1 = gtkVBoxNew();
wizardQ1$packStart(gtkLabelNew("Do you xxxxx ?"));
radioYes = gtkRadioButtonNewWithLabel(NULL, "YES");
radioNo = gtkRadioButtonNewWithLabelFromWidget(radioYes, "NO");
wizardQ1$packStart(radioYes);
wizardQ1$packStart(radioNo);
wizardPage$appendPage(wizardQ1);
wizardPage$setPageTitle(wizardQ1, "Question 1");
wizardPage$setPageComplete(wizardQ1, TRUE);
# Final step
wizardConfirmPage = gtkLabelNew("According to your selections, the suggested template is xxx.
				Press 'Apply' to confirm.");
wizardConfirmPage$modifyFont(pangoFontDescriptionFromString("sans 11"));
wizardConfirmPage$setAlignment(0.2, 0.1);
wizardPage$appendPage(wizardConfirmPage);
wizardPage$setPageTitle(wizardConfirmPage, "Confirmation");
wizardPage$setPageType(wizardConfirmPage, GtkAssistantPageType["confirm"]);
wizardPage$setPageComplete(wizardConfirmPage, TRUE);

# Button box
buttonGroup = ggroup(container = hGroup, horizontal = FALSE);
visible(buttonGroup) = FALSE;
# Buttons
runButton = gbutton("Run", container = buttonGroup);
synButton = gbutton("Syntex check", container = buttonGroup);
anoButton = gbutton("Another", container = buttonGroup);

# Events
# Exit
onExit = function(h, ...) dispose(mainWin);
# Open wizard
onOpenWizard = function(h, ...)
{
	wizardPage$show();
	return(NULL);
}
# Wizard confirmed
onWizardConfirmed = function(widget, param)
{
	widget$hide();
	fileName = "shobbs.Rop";
	filePath = system.file("resources", "Rop", fileName, package = "optimgui");
	if(!length(filePath)) return(NULL);
	openRopFile(filePath);
	return(FALSE);
}
# Open Catalog Page event
onOpenCatalog = function(h, ...)
{
	add(noteBook, catalogPageNew()$catalogPage, label = "Catalog");
	return(NULL);
}
# Default event
onDefaultEvent = function(h, ...) gmessage("Not implemented yet. :(", "Message");
# Open built-in Rop file
onOpenBuiltinRopFile = function(h, ...)
{
	fileName = as.character(svalue(h$action$RopTable));
	filePath = system.file("resources", "Rop", fileName, package = "optimgui");
	if(!length(filePath)) return(NULL);
	openRopFile(filePath);
	return(NULL);
}
# Open Rop file event
onOpenRopFile = function(h, ...)
{
	filePath = gfile(text = "Choose an Rop file...", type = "open",
			         filter = list("Rop files (*.Rop)" = list(patterns = "*.Rop")),
			         handler = function(h,...) return(NULL));
	if(is.na(filePath)) return(NULL);
	openRopFile(filePath);
	return(NULL);
}
# Save Rop file event
onSaveRopFile = function(h, ...)
{
    filePath = gfile(text = "Save an Rop file...", type = "save",
			         filter = list("Rop files (*.Rop)" = list(patterns = "*.Rop")),
			         handler = function(h,...) return(NULL));
	if(is.na(filePath)) return(NULL);
    Encoding(filePath) = "UTF-8";
	writeRop(filePath %+% ".Rop");
	return(NULL);
}
tabsList = list();
# Close Rop file event
onCloseRopFile = function(h, ...)
{
	visible(buttonGroup) = FALSE;
	visible(noteBook) = FALSE;
	while(dispose(noteBook)){}
	currentFile <<- NULL;
	add(noteBook, welcomePageNew()$welcomePage, label = "Welcome");
	visible(noteBook) = TRUE;
    enabled(mSave) = FALSE;
	return(NULL);
}
# Run code event
onRunCode = function(h, ...)
{
	if(!is.null(currentFile))
	{
		z = textConnection("output", open = "w");
		sink(z);
		source(currentFile, echo = FALSE, print.eval = TRUE);
		sink();
		close(z);
		output = paste(output, collapse = "\n");
		outputWidget = tabsList[["Run"]]$rcode;
		svalue(outputWidget) = paste(svalue(outputWidget), output, sep = "\n");
        svalue(noteBook) = which(names(tabsList) == "Run");
	}
}

# Add events
welcomePageAddEvent = function(welcomePage)
{
	addHandlerClicked(welcomePage$wizardLabel, onOpenWizard);
	addHandlerClicked(welcomePage$newLabel, onOpenCatalog);
	addHandlerClicked(welcomePage$openLabel, onOpenRopFile);
}
welcomePage = welcomePageNew();
add(noteBook, welcomePage$welcomePage, label = "Welcome");
catalogPageAddEvent = function(catalogPage)
{
	addHandlerClicked(catalogPage$openRopButton, handler = onOpenBuiltinRopFile,
			action = catalogPage);
}
gSignalConnect(wizardPage, "apply", onWizardConfirmed, NULL);
addHandlerClicked(runButton, onRunCode);


# Menu list
mOpen = gaction(label = "Open", handler = onOpenRopFile, icon = "open");
mSave = gaction(label = "Save", handler = onSaveRopFile, icon = "save");
mClose = gaction(label = "Close", handler = onCloseRopFile, icon = "close");
mExit = gaction(label = "Exit", handler = onExit, icon = "quit");
mCopy = gaction(label = "Copy", handler = onDefaultEvent);
mCut = gaction(label = "Cut", handler = onDefaultEvent);
mTool = gaction(label = "Tool", handler = onDefaultEvent);
mHelp = gaction(label = "optimgui Help", handler = onDefaultEvent);
mAbout = gaction(label = "About", handler = onDefaultEvent);

menuList = list(File = list(open = mOpen, save = mSave, close = mClose,
                            exit = mExit),
		        Edit = list(copy = mCopy, cut = mCut),
		        Tools = list(tool = mTool),
		        Help = list(help = mHelp, about = mAbout));
# Main menu
gmenu(menuList, container = mainWin);
enabled(mSave) = FALSE;
# Show main window
visible(mainWin) = TRUE;
}
