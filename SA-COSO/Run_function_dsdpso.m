%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Implementation of a surrogate-assisted cooperative swarm optimization (SA-COSO) for high-dimensional computationally expensive problems
%%
%%  See the details of SA-COSO in the following paper
%%  C. Sun, Y. Jin, R. Cheng, J. Ding and J. Zeng, Surrogate-assisted Cooperative Swarm Optimization of High-dimensional Expensive Problems,
%%  IEEE Transactions on Evolutionary Computation, 2017
%%
%%  The source code CSO is implemented by Ran Cheng
%%
%%  If you have any questions about the code, please contact:
%%  Chaoli Sun at chaoli.sun.cn@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
for i=1:12
    [gbest_output evaltimes_output final_best max_iteration N_evaluation G_optimum time_output] = MSAMSOA_version2(i);
    if i==1
        save 'd:\test\griewank\gbest_output.txt' gbest_output -ascii;
        save 'd:\test\griewank\evaltimes_output.txt'  evaltimes_output -ascii;
        save 'd:\test\griewank\final_best.txt' final_best -ascii;
        save 'd:\test\griewank\max_iteration.txt' max_iteration -ascii;
        save 'd:\test\griewank\N_evaluation.txt' N_evaluation -ascii;
        save 'd:\test\griewank\G_optimum.txt' G_optimum -ascii;
        save 'd:\test\griewank\time_output.txt' time_output -ascii;
    else if i == 2
            save 'd:\test\test\ackley\gbest_output.txt' gbest_output -ascii;
            save 'd:\test\test\ackley\evaltimes_output.txt'  evaltimes_output -ascii;
            save 'd:\test\test\ackley\final_best.txt' final_best -ascii;
            save 'd:\test\test\ackley\max_iteration.txt' max_iteration -ascii;
            save 'd:\test\test\ackley\N_evaluation.txt' N_evaluation -ascii;
            save 'd:\test\test\ackley\G_optimum.txt' G_optimum -ascii;
            save 'd:\test\test\ackley\time_output.txt' time_output -ascii;
        else if i == 3
                save 'd:\test\test\rosenbrock\gbest_output.txt' gbest_output -ascii;
                save 'd:\test\rosenbrock\evaltimes_output.txt'  evaltimes_output -ascii;
                save 'd:\test\rosenbrock\final_best.txt' final_best -ascii;
                save 'd:\test\rosenbrock\max_iteration.txt' max_iteration -ascii;
                save 'd:\test\rosenbrock\N_evaluation.txt' N_evaluation -ascii;
                save 'd:\test\rosenbrock\G_optimum.txt' G_optimum -ascii;
                save 'd:\test\rosenbrock\time_output.txt' time_output -ascii;
            else if i == 4
                    save 'd:\test\ellipsoid\gbest_output.txt' gbest_output -ascii;
                    save 'd:\test\ellipsoid\evaltimes_output.txt'  evaltimes_output -ascii;
                    save 'd:\test\ellipsoid\final_best.txt' final_best -ascii;
                    save 'd:\test\ellipsoid\max_iteration.txt' max_iteration -ascii;
                    save 'd:\test\ellipsoid\N_evaluation.txt' N_evaluation -ascii;
                    save 'd:\test\ellipsoid\G_optimum.txt' G_optimum -ascii;
                    save 'd:\test\ellipsoid\time_output.txt' time_output -ascii;
                else if i == 5
                        save 'd:\test\griewank100\gbest_output.txt' gbest_output -ascii;
                        save 'd:\test\griewank100\evaltimes_output.txt'  evaltimes_output -ascii;
                        save 'd:\test\griewank100\final_best.txt' final_best -ascii;
                        save 'd:\test\griewank100\max_iteration.txt' max_iteration -ascii;
                        save 'd:\test\griewank100\N_evaluation.txt' N_evaluation -ascii;
                        save 'd:\test\griewank100\G_optimum.txt' G_optimum -ascii;
                        save 'd:\test\griewank100\time_output.txt' time_output -ascii;
                    else if i == 6
                            save 'd:\test\ackley100\gbest_output.txt' gbest_output -ascii;
                            save 'd:\test\ackley100\evaltimes_output.txt'  evaltimes_output -ascii;
                            save 'd:\test\ackley100\final_best.txt' final_best -ascii;
                            save 'd:\test\ackley100\max_iteration.txt' max_iteration -ascii;
                            save 'd:\test\ackley100\N_evaluation.txt' N_evaluation -ascii;
                            save 'd:\test\ackley100\G_optimum.txt' G_optimum -ascii;
                            save 'd:\test\ackley100\time_output.txt' time_output -ascii;
                        else if i == 7
                                save 'd:\test\FESPSO_RBF1\dsdpso\rosenbrock100\gbest_output.txt' gbest_output -ascii;
                                save 'd:\test\FESPSO_RBF1\dsdpso\rosenbrock100\evaltimes_output.txt'  evaltimes_output -ascii;
                                save 'd:\test\FESPSO_RBF1\dsdpso\rosenbrock100\final_best.txt' final_best -ascii;
                                save 'd:\test\FESPSO_RBF1\dsdpso\rosenbrock100\max_iteration.txt' max_iteration -ascii;
                                save 'd:\test\FESPSO_RBF1\dsdpso\rosenbrock100\N_evaluation.txt' N_evaluation -ascii;
                                save 'd:\test\FESPSO_RBF1\dsdpso\rosenbrock100\G_optimum.txt' G_optimum -ascii;
                                save 'd:\test\FESPSO_RBF1\dsdpso\rosenbrock100\time_output.txt' time_output -ascii;
                            else if i == 8
                                    save 'd:\test\FESPSO_RBF1\dsdpso\ellipsoid100\gbest_output.txt' gbest_output -ascii;
                                    save 'd:\test\FESPSO_RBF1\dsdpso\ellipsoid100\evaltimes_output.txt'  evaltimes_output -ascii;
                                    save 'd:\test\FESPSO_RBF1\dsdpso\ellipsoid100\final_best.txt' final_best -ascii;
                                    save 'd:\test\FESPSO_RBF1\dsdpso\ellipsoid100\max_iteration.txt' max_iteration -ascii;
                                    save 'd:\test\FESPSO_RBF1\dsdpso\ellipsoid100\N_evaluation.txt' N_evaluation -ascii;
                                    save 'd:\test\FESPSO_RBF1\dsdpso\ellipsoid100\G_optimum.txt' G_optimum -ascii;
                                    save 'd:\test\FESPSO_RBF1\dsdpso\ellipsoid100\time_output.txt' time_output -ascii;
                                else if i == 9
                                        save 'd:\test\f10\gbest_output.txt' gbest_output -ascii;
                                        save 'd:\test\f10\evaltimes_output.txt'  evaltimes_output -ascii;
                                        save 'd:\test\f10\final_best.txt' final_best -ascii;
                                        save 'd:\test\f10\max_iteration.txt' max_iteration -ascii;
                                        save 'd:\test\f10\N_evaluation.txt' N_evaluation -ascii;
                                        save 'd:\test\f10\G_optimum.txt' G_optimum -ascii;
                                        save 'd:\test\f10\time_output.txt' time_output -ascii;
                                    else if i == 10
                                            save 'd:\test\f10_100\gbest_output.txt' gbest_output -ascii;
                                            save 'd:\test\f10_100\evaltimes_output.txt'  evaltimes_output -ascii;
                                            save 'd:\test\f10_100\final_best.txt' final_best -ascii;
                                            save 'd:\test\f10_100\max_iteration.txt' max_iteration -ascii;
                                            save 'd:\test\f10_100\N_evaluation.txt' N_evaluation -ascii;
                                            save 'd:\test\f10_100\G_optimum.txt' G_optimum -ascii;
                                            save 'd:\test\f10_100\time_output.txt' time_output -ascii;
                                        else if i == 11
                                                save 'd:\test\f19\gbest_output.txt' gbest_output -ascii;
                                                save 'd:\test\f19\evaltimes_output.txt'  evaltimes_output -ascii;
                                                save 'd:\test\f19\final_best.txt' final_best -ascii;
                                                save 'd:\test\f19\max_iteration.txt' max_iteration -ascii;
                                                save 'd:\test\f19\N_evaluation.txt' N_evaluation -ascii;
                                                save 'd:\test\f19\G_optimum.txt' G_optimum -ascii;
                                                save 'd:\test\f19\time_output.txt' time_output -ascii;
                                            else if i == 12
                                                    save 'd:\test\f19_100\gbest_output.txt' gbest_output -ascii;
                                                    save 'd:\test\f19_100\evaltimes_output.txt'  evaltimes_output -ascii;
                                                    save 'd:\test\f19_100\final_best.txt' final_best -ascii;
                                                    save 'd:\test\f19_100\max_iteration.txt' max_iteration -ascii;
                                                    save 'd:\test\f19_100\N_evaluation.txt' N_evaluation -ascii;
                                                    save 'd:\test\f19_100\G_optimum.txt' G_optimum -ascii;
                                                    save 'd:\test\f19_100\time_output.txt' time_output -ascii;
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
