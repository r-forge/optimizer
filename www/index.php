
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='http://r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<!-- R-Forge Logo -->
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="/"><img src="<?php echo $themeroot; ?>/images/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<?php if ($handle=fopen('http://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')){
$contents = '';
while (!feof($handle)) {
	$contents .= fread($handle, 8192);
}
fclose($handle);
echo $contents; } ?>

<!-- end of project description -->

<p> The <strong>project summary page</strong> you can find <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. </p>

<p><p> <strong>BB</strong> -  Solving and optimizing large-scale nonlinear systems
<p> Barzilai-Borwein Spectral Methods for solving nonlinear system of 
 equations, and for optimizing nonlinear objective functions subject to simple 
 constraints.

<p><p><strong>GPArotation</strong> - GPA Factor Rotation
<p><p>Gradient Projection Algorithm Rotation for Factor Analysis. 

<p><p><strong>numDeriv</strong> - Accurate Numerical Derivatives
<p>Accurate Numerical Derivatives. See ?numDeriv.Intro for more details.

<p><p><strong>optimx</strong> - A wrapper to integrate optim() and related optimization functionality
<p>The package allows users to conduct optimization of smooth, continuous functions 
using the methods in optim (Nelder-Mead, BFGS, L-BFGS-B, CG, and SANN) as well as
those from packages nlm, nlminb, and ucminf. A number of convenience features are
available. At 2009-06-24, this is still a beta level package.</p>

<p><p><strong>funcheck</strong> - Checking functions that are to be optimized
<p>The package consists of functions funcheck() to verify files containing the function and
various derivative computations in a standardized naming convention. Function funtest carries
out the same tests on explicitly named objective and derivative codes.</p>

</body>
</html>
