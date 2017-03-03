process foo {
  tag { ${code.replaceAll(/''/, fold)} }

  input:
  val code from 'alpha', 'gamma', 'omega'
  val fold from 'dalpha', 'dgamma', 'domega'
  
  """
  echo $code
  """
}
