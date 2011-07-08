optimgui = function(){

# Main Window
mainWin = gwindow("optimgui", visible = FALSE);
size(mainWin) = c(800, 600);
# Horizontal layout box
hGroup = ggroup(container = mainWin, expand = TRUE);
# The editor
editor <<- EditorNew(container = hGroup);

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
buttonGroup <<- ggroup(container = hGroup, horizontal = FALSE);
visible(buttonGroup) = FALSE;
# Buttons
runButton = gbutton("Run", container = buttonGroup);
synButton = gbutton("Syntex check", container = buttonGroup);
anoButton = gbutton("Another", container = buttonGroup);

# Add events
editor = showWelcomePage.Editor(editor);
gSignalConnect(wizardPage, "apply", onWizardConfirmed, NULL);
addHandlerClicked(runButton, onRunCode);


# Menu list
mOpen = gaction(label = "Open", handler = onOpenRopFile, icon = "open");
mSave <<- gaction(label = "Save", handler = onSaveRopFile, icon = "save");
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
