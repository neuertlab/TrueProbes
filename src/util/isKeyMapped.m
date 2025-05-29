function out = isKeyMapped(key,Map)
try
   out = Map(key); 
catch
   out = 0; 
end
end