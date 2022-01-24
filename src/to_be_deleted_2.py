

class A:
    def __init__(self) -> None:
        self.x = 1
        self.y = 2
    
    def test(self, data):
        print(id(data))


class B:
    def __init__(self, other) -> None:
        self.x = other.x
        self.y = other.y
        other.x = 10

foo = A()

print(id(foo.x))
foo.test(foo.x)

bar = B(foo)

print(foo.x)
print(foo.__dict__)
print(bar.__dict__)