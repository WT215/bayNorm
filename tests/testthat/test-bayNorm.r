set.seed(12300)
data("EXAMPLE_DATA_list")

# bay_array_N<-array(rpois(1000*50*2,17),dim=c(1000,50,2))
# bay_array_C<-array(rpois(1000*50*2,58),dim=c(1000,50,2))
# \dontrun{
# noisy_output<-NOISY_FUN(bay_array_N,bay_array_C)

test_that("bayNorm runs properly",
          {
set.seed(1258484)
q1<- bayNorm(Data=EXAMPLE_DATA_list$inputdata[seq(1,10),],BETA_vec = EXAMPLE_DATA_list$inputbeta,parallel = F)
set.seed(1258484)
q2<- bayNorm(Data=EXAMPLE_DATA_list$inputdata[seq(1,10),],BETA_vec = EXAMPLE_DATA_list$inputbeta,parallel = F)
expect_identical(q1, q2)
}
)


test_that("bayNorm_sup runs properly",
          {
            set.seed(1258484)
            q1<- bayNorm(Data=EXAMPLE_DATA_list$inputdata[seq(1,10),],BETA_vec = EXAMPLE_DATA_list$inputbeta,parallel = F)

            set.seed(1258484)
            q2<- bayNorm_sup(Data=EXAMPLE_DATA_list$inputdata[seq(1,10),],BETA_vec = EXAMPLE_DATA_list$inputbeta,PRIORS=q1$PRIORS,parallel = F)
            expect_identical(q1$PRIORS, q2$PRIORS)
            }
)

