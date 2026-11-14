# Runtime threading policy. Julia's thread pool is fixed at startup, but every
# PIC loop consults this flag and can therefore run serially on another machine.
const PIC_THREADING_ENABLED = Ref(false)

default_machine_threading_policy() = (
    automatic=true,
    enable_on_matching_machine=true,
    required_hostname="Scarlett-Z.local",
    required_cpu_model="Apple M5",
    required_hardware_cores=10,
    threads=10,
    manual_enabled=false,
)

macro pic_threads(loop)
    threaded_loop = :(Threads.@threads $loop)
    return esc(quote
        if PIC_THREADING_ENABLED[]
            $threaded_loop
        else
            $loop
        end
    end)
end

function configure_machine_threading!(policy)
    cpu_models = unique(info.model for info in Sys.cpu_info())
    detected_model = isempty(cpu_models) ? "unknown" : join(cpu_models, ", ")
    detected_cores = length(Sys.cpu_info())
    detected_hostname = gethostname()
    matches = Sys.isapple() && occursin("apple", lowercase(Sys.MACHINE)) &&
              detected_hostname == policy.required_hostname &&
              detected_model == policy.required_cpu_model &&
              detected_cores == policy.required_hardware_cores

    PIC_THREADING_ENABLED[] = policy.automatic ?
        (policy.enable_on_matching_machine && matches) : policy.manual_enabled

    if PIC_THREADING_ENABLED[]
        println("Machine check passed: $detected_hostname, $detected_model, " *
                "$detected_cores hardware cores. PIC multithreading enabled " *
                "with $(Threads.nthreads()) Julia thread(s).")
        Threads.nthreads() != policy.threads &&
            println("  Warning: configured $(policy.threads) Julia threads, but " *
                    "this process has $(Threads.nthreads()). Use scripts/run_pic.jl " *
                    "for automatic relaunch.")
    else
        println("Machine check did not match the configured MacBook; " *
                "PIC kernels will run serially. Detected: $detected_hostname, " *
                "$detected_model, $detected_cores hardware cores, $(Sys.KERNEL).")
    end
    return PIC_THREADING_ENABLED[]
end
