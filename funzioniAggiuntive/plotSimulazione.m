function [plot_T , plot_Q , plot_U]  = plotSimulazione(tempo, storia_x  , storia_u , x_ref)
    tempo = tempo/60; %[min]

    figure
    
    sgtitle("Evoluzioni degli stati")
    
    plot_T = subplot(2 , 1, 1);
    plot( tempo, storia_x(1:3 , :))
    yline(x_ref(1))
    legend(["T1" , "T2" , "T3" ,"Obbiettivo"])
    ylabel("Temperatura $[^{\circ}C]$" , Interpreter="latex");
    xlabel("Tempo $[min]$" , Interpreter="latex");
    title("Temperatura")
    
    plot_Q = subplot(2 , 1, 2);
    plot(tempo, storia_x(4:end , :))
    yline(x_ref(4))
    ylim([0 , 150])
    legend(["Q1"  "Q2"  "Q3" "Obbiettivo"])
    ylabel("Potenza $[W]$" , Interpreter="latex");
    xlabel("Tempo $[min]$" , Interpreter="latex");
    title("Potenza dei termosifoni")
    
    figure
    plot_U = plot(tempo, storia_u);
    title("Azioni di controllo")
    ylabel("Potenza $[W]$" , Interpreter="latex");
    xlabel("Tempo $[min]$" , Interpreter="latex");
    legend(["Q1"  "Q2"  "Q3"])

end