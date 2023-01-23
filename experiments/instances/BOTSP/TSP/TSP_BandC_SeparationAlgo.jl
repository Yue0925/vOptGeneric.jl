using JuMP
using CPLEX
using CPUTime

function ViolatedMengerCut_IntegerSeparation(G,xsep)

  #list<int>::const_iterator it;
  #list<C_link *>::iterator ita;
 

  start=rand(1:G.nb_points);
  
  W =find_cycle_in_integer_x(xsep, start)

  if size(W)!=G.nb_points
      
      con = @build_constraint(sum(x[i,j] for i ∈ W for j ∈ 1:G.nb_points if j ∉ W) >= 2)
      
      return true, con
  else
      return false , @build_constraint(0>=0)
  end


end
