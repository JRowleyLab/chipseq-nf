process sayHello {
    input:
        val cheers
    output:
        stdout

    """
    echo $cheers
    """
}