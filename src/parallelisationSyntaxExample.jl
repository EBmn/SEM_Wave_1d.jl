

function main() 

    # this code will not run
    
    omegas = [1, 2, 3]    # some parameter values for which I will run my simulations
    
    n = length(omegas)   # number of tasks that will be completed
    ch = Channel{Int}(n) # make channel that will keep track of tasks
                         
    for j in eachindex(omegas)      # fill channel with work, in this case indices of the list omegas
        put!(ch, j)
    end
    
    close(ch)            # all tasks have been added to the channel, so we close it
    
    

    completed = falses(n) # one solution if you want to keep track of which tasks have been completed already
    
    tasks = Vector{Task}(undef, Threads.nthreads()) # not entirely sure whether this is necessary, or whether it is an artefact of a previous version of the code. In any case, this works for me.

    for t in 1:length(tasks)

        tasks[t] = Threads.@spawn begin # spawn a thread... again I don't think I actually need to store the threads like this anymore.
        
    # I'm guessing something like this also works:
    #=
    # substitute lines 24-26 for this
    nT = Threads.nthreads() # number of threads available to Julia right now
    for t = 1:nT 
        Threads.@spawn begin
    =#
            
            try # ... and put it to work. The try/catch structure is not necessary, but nice for debugging if something crashes.
                for j in ch # the thread takes work from the channel ch (in this case we have stored indices of the list of parameter values omegas)

                    
                    tid = threadid()
                    omega = omegas[j]
                    
                    #println("thread $tid considering value $omega")
                    data = simulate(omega) # or whatever, do the actual work here!

                    # Now do whatever we want with the data once we have it. Perhaps store it somewhere, but be careful about deadlocks, race conditions, or similar.

                    completed[j] = true # keep track of which tasks have been completed
                        
                end
                    
            
            catch e

                # if something goes wrong with a thread, we print the error message here
                
                println("Worker thread $(threadid()) crashed!")
                showerror(stdout, e)
                println()
                rethrow(e)

            end
            
        end

    end


    wait.(tasks) # wait until all the tasks are done


    # do whatever we want with our output here

end
    
