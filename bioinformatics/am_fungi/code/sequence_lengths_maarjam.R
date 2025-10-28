

data<-readDNAStringSet("data/ref_seqs/maarjam_SSU_NS31_2021.fasta")

# Extract sequence lengths
seq_lengths <- width(data)

# Create a data frame for plotting
df <- data.frame(length = seq_lengths)

# Plot density using ggplot2
ggplot(df, aes(x = length)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  labs(title = "Density Plot of Sequence Lengths",
       x = "Sequence Length",
       y = "Density") +
  theme_minimal()

# Get the sequences less than 440 in length
short_seqs <- data[seq_lengths < 440]

# Make a tibble of the names and lengths of these short sequences
short_seqs_df <- tibble(
  name = names(short_seqs),
  length = width(short_seqs)
)

# MAke a tibble of the names and lengths of the sequences greater than 440
long_seqs_df <- tibble(
  name = names(data[seq_lengths >= 440]),
  length = width(data[seq_lengths >= 440])
)

short_seqs_df %>%
  filter(grepl(" Glomus ", name))

long_seqs_df %>%
  filter(grepl(" Glomus ", name))
