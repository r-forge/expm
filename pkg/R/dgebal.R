dgebal <- function(A, job = c("B","P","S"))
{
    job <- match.arg(job)
    n <- dim(A)[1]
    .Call("R_dgebal", A,  job)
}
