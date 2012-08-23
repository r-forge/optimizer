
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

<p><p> Packages for optimization and solution of nonlinear systems of equations include
<ul>
<li>- large scale optimization with box constraints</li>
<li>- exterior point method for general constraint optimization</li>
<li>- sequential quadratic programming with box constraints</li>
<li>- conjugate gradient optimization with box constraints</li>
<li>- fixed point methods</li>
<li>- a Levenberg-Marquardt approach to nonlinear least squares</li>
<li>- large scale equation solving via Barzilai-Borwein Spectral Methods for solving nonlinear system of 
 equations</li>
<li>- multi-objective optimization for Nash equilibrium</li>
</ul>


<p><p> Packages for utils functions include
<ul>
<li>- KKT condition testing</li>
<li>- numerical derivative routine</li>
<li>- a common interface for optimization methods</li>
</ul> 
 

<h3>Extra documentation</h3>

<p>There is an attempt to build an evolving "handbook" or "cheatsheet" 
for nonlinear least squares computations for <strong>R</strong>. The 
source file is of type .Rnw and can be processed with the <strong>knitR</strong>
tools. This is <a href="http://optimizer.r-forge.r-project.org/nlshb.Rnw">nlshb.Rnw</a>, with the latest
output file uploaded as <a href="http://optimizer.r-forge.r-project.org/nlshb.pdf">nlshb.pdf</a>


</body>
</html>
