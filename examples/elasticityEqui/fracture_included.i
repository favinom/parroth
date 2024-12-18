[Mesh]
 [./cmg]
    type = CartesianMeshGenerator
    dim = 2
    dx = '0.4995 1e-3 0.4995'
    dy = '0.25 0.5 0.25'
    ix = '10 2 10'
    iy = '5 10 5'
    subdomain_id = '0 0 0 0 1 0 0 0 0'
  [../]
[]

[Variables]
  [disp_x] []
  [disp_y] []
[]

[Kernels]
 [lex] type = LinearElasticity variable = disp_x disp_x = disp_x disp_y = disp_y component = 0 []
 [ley] type = LinearElasticity variable = disp_y disp_x = disp_x disp_y = disp_y component = 1 []
[]

[Materials]
 [./ElasticityMaterialProperties0] type = ElasticityMaterialProperties mu  = 1.0   lambda = 1.0 disp_x = disp_x disp_y = disp_y block = 0 [../]
 [./ElasticityMaterialProperties1] type = ElasticityMaterialProperties mu  = 0.001 lambda = 0.01 disp_x = disp_x disp_y = disp_y block = 1 [../]
[]


[BCs]
 [left_x] type = DirichletBC value = 0 variable = disp_x boundary = left []
 [left_y] type = DirichletBC value = 0 variable = disp_y boundary = left []
 [right_x] type = DirichletBC value = 0.01 variable = disp_x boundary = right []
 [right_y] type = DirichletBC value = 0.01 variable = disp_y boundary = right []
[]


[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  line_search = none
#  nl_rel_tol = 1e-10
#  nl_abs_tol = 1e-10

# petsc_options = '-ksp_view_pmat'
 petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
 petsc_options_value='  preonly   lu       NONZERO               mumps '

# [Quadrature] type = TRAP []

[]

[Outputs]
  file_base = equi_fracture_included
  print_linear_residuals = false
  exodus = true
#  xda = true
[]

