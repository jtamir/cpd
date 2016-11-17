
cpdsrcs := $(wildcard $(srcdir)/cpd/*.c)
cpdobjs := $(cpdsrcs:.c=.o)

.INTERMEDIATE: $(cpdobjs)

lib/libcpd.a: libcpd.a($(cpdobjs))



