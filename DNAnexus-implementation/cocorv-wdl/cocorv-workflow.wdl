version 1.0 

workflow hello_world { 
    input { 
        String name 
    }

    call write_greeting { 
        input: greet_name = name 
    }
}

task write_greeting { 
    input {
        String greet_name 
    }

    command { 
        echo 'Hello, ${greet_name}!' > out.txt
    }

    output {
        File outfile = "out.txt" 
    }
}