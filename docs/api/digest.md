# Digestion and filtering

`hicstuff.digest` cuts a genome into fragments (`write_frag_info`,
`get_restriction_table`) and assigns the mates of each pair to them
(`attribute_fragments`). `hicstuff.filter` finds the uncuts and loops of a
library (`get_thresholds`) and removes them (`filter_events`).

See [Digesting a genome](../guides/digest.md) and
[Filtering 3C events](../guides/filter.md).

::: hicstuff.digest

::: hicstuff.filter
