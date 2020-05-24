function Cpm=fCpm(Tm)
if Tm<318
    Cpm=2000;
else
    if 318<Tm && Tm<320
        Cpm=1e5;
    else
        Cpm=2010;
    end
end

end
