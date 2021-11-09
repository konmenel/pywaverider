

class A:
    def __init__(self) -> None:
        self.x = 1
        self.y = 2


class B:
    def __init__(self, other) -> None:
        self.x = other.x
        self.y = other.y
        other.x = 10

foo = A()
bar = B(foo)

print(foo.x)
print(foo.__dict__)
print(bar.__dict__)